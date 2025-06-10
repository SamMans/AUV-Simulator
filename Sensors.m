         % Author: Samer A. Mohamed (https://github.com/SamMans) %
classdef Sensors
    properties(SetAccess = private)
    % Pressure sensor parameters
        dT_p %--> Sensor refresh sampling time
        TRANS_p %--> Pressure sensor mounting shift from COM
        Stdev_p %--> Pressure sensor reading standard deviation
        
    % IMU parameters
        dT_imu %--> Sensor refresh sampling time
        Stdev_imu %--> IMU standard deviation
        
    % Stereo camera parameters
        dT_ster %--> Stereo camera most expected refresh rate
        Stdev_ster %--> Depth data standard deviation
        TRANS_ster %--> Camera transformation from Camera frame to COM frame
        Shfov %--> Camera horizontal field of view
        Svfov %--> Camera vertical field of view
        Srange %--> Camera depth range 
        LM %--> Global positions of landmarks in environment (x,y,z)
    end
    methods
        function obj = Sensors(models, Pshift, Sshift, Landmarks) %--> Constructor
            % Initialize pressure sensor parameters
            obj.TRANS_p = Pshift;
            load(pwd + "\Libraries\Sensors\Depth sensors\" + models(1) + ".mat");
            obj.Stdev_p = Pprec / 300; obj.dT_p = 1 / Fp;
            
            % Initialize IMU parameters
            load(pwd + "\Libraries\Sensors\IMUs\" + models(2) + ".mat");
            obj.Stdev_imu = imu_dev; obj.dT_imu = 1 / Fimu;
            
            % Initialize stereo camera parameters
            obj.LM = Landmarks;
            load(pwd + "\Libraries\Sensors\Stereo cameras\" + models(3) + ".mat");
            obj.TRANS_ster = inv([1 0 0 Sshift(1);         
               0 1 0 Sshift(2);
               0 0 1 Sshift(3);
               0 0 0 1]*...
                [cos(pi/2) 0 sin(pi/2) 0;
                   0 1 0 0;
               -sin(pi/2) 0 cos(pi/2) 0;
                   0 0 0 1]*...
              [cos(pi/2) -sin(pi/2) 0 0;
                sin(pi/2) cos(pi/2) 0 0;
                     0 0 1 0;
                     0 0 0 1]); %--> Camera matrix (COM w.r.t Camera)
            obj.dT_ster = 1 / Fs; obj.Stdev_ster = prec_s / 300;
            obj.Shfov = Hfov; obj.Svfov = Vfov; obj.Srange = range_s;
        end
        function J = Jac(obj, n)
            % Jacobian of the transform (COM frame w.r.t global frame)
            J = [cos(n(6))*cos(n(5)) -sin(n(6))*cos(n(4))+cos(n(6))*sin(n(5))*sin(n(4)) sin(n(6))*sin(n(4))+cos(n(6))*cos(n(4))*sin(n(5)) n(1);
            sin(n(6))*cos(n(5)) cos(n(6))*cos(n(5))+sin(n(4))*sin(n(5))*sin(n(6)) -cos(n(6))*sin(n(4))+sin(n(5))*sin(n(6))*cos(n(4)) n(2);
            -sin(n(5)) cos(n(5))*sin(n(4)) cos(n(5))*cos(n(4)) n(3);
            0 0 0 1];
        end
        function [P_read, imu_read, ster_read] = Refresh(obj, n, Rt)
            % Pressure sensor refresh (must follow standard refresh rate)
            if(rem(round(Rt,1), obj.dT_p) == 0)
                % Apply homogenous transformation
                R = [1 0 0 0;
                    0 cos(n(4)) sin(n(4)) 0;
                    0 -sin(n(4)) cos(n(4)) 0;
                    0 0 0 1]; %--> Roll rotation
                P = [cos(n(5)) 0 -sin(n(5)) 0;
                    0 1 0 0;
                    sin(n(5)) 0 cos(n(5)) 0;
                    0 0 0 1]; %--> Pitch rotation
                Y = [cos(n(6)) sin(n(6)) 0 0;
                    -sin(n(6)) cos(n(6)) 0 0;
                     0 0 1 0;
                     0 0 0 1]; %--> Yaw rotation
                Tr = Y * P * R; %--> Full transform
                Shift = (Tr * [obj.TRANS_p, 1].').';
                
                % Compute final sensor reading
                P_read = n(3) + Shift(3) + normrnd(0, obj.Stdev_p); %--> Add noise & shift
            else
                P_read = []; %--> No reading currently
            end
            
            % IMU sensor refresh (must follow standard refresh rate)
            if(rem(round(Rt,1), obj.dT_imu) == 0)
                imu_read = [n(4), n(5), n(6)] + normrnd(0, obj.Stdev_imu); %--> Angle readings + noise
            else
                imu_read = []; %--> No reading currently
            end
            
            % Stereo camera refresh (follows a periodic randomization)
            if(rem(round(Rt,1), obj.dT_ster) == 0)
                % Create the transform "Global frame w.r.t camera frame"
                T_cam = obj.TRANS_ster * inv(Jac(obj, n));
                
                % Initialize an empty message
                ster_read = nan * ones(size(obj.LM, 1), 4);
                
                % Loop through landmarks to see which one lies in f.o.v
                for i = 1 : size(obj.LM, 1)
                    % Transform locations of landmarks to camera frame
                    loc = (T_cam * [obj.LM(i, :), 1].').';
                    
                    % Check if the landmarks lie within the cone of the camera
                    % If so, camera captures the landmark
                    if(loc(3) <= obj.Srange(2) && loc(3) >= obj.Srange(1) && ...
                            abs(loc(1)) <= loc(3) * tan(obj.Shfov / 2) && ...
                            abs(loc(2)) <= loc(3) * tan(obj.Svfov / 2))
                        % Add sensor noise
                        loc = loc(1:3) + normrnd(0, obj.Stdev_ster);
                        
                        % Append readings to perception
                        ster_read(i, :) = [i-1, loc];
                    end
                end
                
                % Eliminate nan rows (rows without readings)
                ster_read(any(isnan(ster_read), 2), :) = [];
            else
                ster_read = []; %--> No reading currently
            end
        end
    end
end