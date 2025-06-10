         % Author: Samer A. Mohamed (https://github.com/SamMans) %
classdef AUV
    properties(SetAccess = private)
    % Static properties
        m %--> AUV mass
        V %--> AUV volume
        I0 %--> AUV moments of inertia matrix
        M %--> Inertia tensor
        COB %--> Center of buoyancy vector (referenced to COM)
        CD %--> Drag coefficients matrix
        MA %--> Added mass matrix
        CF %--> Coefficient of lift matrix
        proj %--> Thruster projection matrix
        W %--> Weight force
        B %--> Buoyancy force
        Tht %--> Thruster mounting angle (w.r.t surge axis)
        
    % Dynamic properties
        n %--> Global pose vector
        ndot %--> Global pose rate of change
        v %--> Local velocity
        vdot %--> Local acceleration
        tt %--> Thruster current thrust values
        time %--> Time in seconds
        iter %--> Simulation iteration
        dT_min %--> Minimum time step size
        
    % AUV thrusters
        Th
        
    % AUV sensors (pressure sensor, IMU, stereo camera)
        Sens
    end
    methods
        function obj = AUV(model, LM, ni, vi, vdoti) %---> Constructor
            % Read AUV parameters from .mat file + Identify usual constants
            if ~exist(pwd + "\Libraries\AUVs\" + model, 'dir')
               % Throw error if user chooses an unavailable model
               error(model + " doesn't exist, check the available models on SamSim");
            end
            load(pwd + "\Libraries\AUVs\" + model + "\physical.mat");
            g = 9.81; raw = 1000;
            
            % Create a T200 object (8 Thrusters)
            obj.Th = [Thruster(1, Th_mod(1), Th_mod(2)),Thruster(2, Th_mod(1), Th_mod(2)), ...
                Thruster(3, Th_mod(1), Th_mod(2)), Thruster(4, Th_mod(1), Th_mod(2)), ...
                Thruster(5, Th_mod(1), Th_mod(2)), Thruster(6, Th_mod(1), Th_mod(2)), ...
                Thruster(7, Th_mod(1), Th_mod(2)), Thruster(8, Th_mod(1), Th_mod(2))]; 
            
            % Create a sensors object
            obj.Sens = Sensors(Sens_mod, Pshift, Sshift, LM);
            
            % Calculate Hydrostatic parameters
            obj.m = m_AUV; obj.V = V_AUV; obj.W = obj.m*g; 
            obj.B = raw*obj.V*g; obj.I0 = I0_AUV; obj.COB = COB_AUV;
            obj.M = [diag([obj.m,obj.m,obj.m]) zeros(3,3);
               zeros(3,3)       obj.I0 ];  %--> Inertia Matrix
           
            % Calculate Hydrodynamic parameters
            obj.MA = MA_AUV; %--> Added mass inertia matrix
            obj.CD = CD_AUV; %--> Drag matrix
            obj.CF = CF_AUV; %--> Lift matrix
            
            % Calculate Thrust forces
            Tht = Tht_AUV; R_ind = R_ind_AUV; Rx = Rx_AUV;
            Ry = Ry_AUV; R_Hor = R_Hor_AUV;
            obj.proj = [cos(Tht),cos(Tht),cos(Tht),cos(Tht),0,0,0,0;
                        -sin(Tht),sin(Tht),-sin(Tht),sin(Tht),0,0,0,0;
                        0,0,0,0,1,1,1,1;
                        sin(Tht)*R_ind,-sin(Tht)*R_ind,sin(Tht)*R_ind,-sin(Tht)*R_ind,Ry,-Ry,-Ry,Ry;
                        cos(Tht)*R_ind,cos(Tht)*R_ind,cos(Tht)*R_ind,cos(Tht)*R_ind,-Rx,-Rx,Rx,Rx;
                        -R_Hor,R_Hor,R_Hor,-R_Hor,0,0,0,0]; %--> Projection matrix
                    
             % Initial conditions
             obj.n = ni; %--> global states
             obj.v = vi; %--> local states
             obj.vdot = vdoti; %--> local states rate of change
             obj.ndot = (Jacobian(obj)*vi.').'; %--> global states rate of change
             obj.time = 0; %--> AUV clock
             obj.iter = 0; %--> State update iterations
             obj.dT_min = 10^-4; %--> Minumum step size 
        end
        function N = getn(obj)
            N = obj.n; %--> get AUV current global states
        end
        function V = getv(obj)
            V = obj.v; %--> get AUV current local states
        end
        function Vdot = getvdot(obj)
            Vdot = obj.vdot; %--> get AUV current local states r.o.c
        end
        function obj = ESC(obj,PWM)
            % Write PWM values to thruster Escs & Update thrust values
            obj.Th(1) = P2T(obj.Th(1),PWM(1));obj.Th(2) = P2T(obj.Th(2),PWM(2));
            obj.Th(3) = P2T(obj.Th(3),PWM(3));obj.Th(4) = P2T(obj.Th(4),PWM(4));
            obj.Th(5) = P2T(obj.Th(5),PWM(5));obj.Th(6) = P2T(obj.Th(6),PWM(6));
            obj.Th(7) = P2T(obj.Th(7),PWM(7));obj.Th(8) = P2T(obj.Th(8),PWM(8));
        end
        function T = getTh(obj)
            T = [getT(obj.Th(1)),getT(obj.Th(2)),getT(obj.Th(3))...
                ,getT(obj.Th(4)),getT(obj.Th(5)),getT(obj.Th(6)),...
                getT(obj.Th(7)),getT(obj.Th(8))]; %--> get AUV current thruster forces 
        end
        function t = gett(obj)
            t = obj.time; %--> Return AUV clock
        end
        function i = geti(obj)
            i = obj.iter; %--> Return AUV simulation iteration
        end
        function dt_min = getdT_min(obj)
            dt_min = obj.dT_min; %--> Return simulation minimum step size
        end
        function J = Jacobian(obj)
            % Return jacobian relating global state rate of change and local velocity
            J = [cos(obj.n(6))*cos(obj.n(5)) -sin(obj.n(6))*cos(obj.n(4))+cos(obj.n(6))*sin(obj.n(5))*sin(obj.n(4))...
                sin(obj.n(6))*sin(obj.n(4))+cos(obj.n(6))*cos(obj.n(4))*sin(obj.n(5)) 0 0 0;
                sin(obj.n(6))*cos(obj.n(5)) cos(obj.n(6))*cos(obj.n(5))+sin(obj.n(4))*sin(obj.n(5))*sin(obj.n(6)) ...
                -cos(obj.n(6))*sin(obj.n(4))+sin(obj.n(5))*sin(obj.n(6))*cos(obj.n(4)) 0 0 0;
                -sin(obj.n(5)) cos(obj.n(5))*sin(obj.n(4)) cos(obj.n(5))*cos(obj.n(4)) 0 0 0;
                0 0 0 1 sin(obj.n(4))*tan(obj.n(5)) cos(obj.n(4))*tan(obj.n(5))
                0 0 0 0 cos(obj.n(4)) -sin(obj.n(4));
                0 0 0 0 sin(obj.n(4))/cos(obj.n(5)) cos(obj.n(4))/cos(obj.n(5))];
        end
        function [obj] = Newton(obj,PWM)
            % Determine current drag, lift and added mass matrices 
            dir = zeros(6,6); 
            for i = 1:6
                dir(i,:) = obj.v >= 0; %--> Motion direction (+ve or -ve)
            end
            cd = dir.*obj.CD(1:6,:) + (1 - dir).*obj.CD(7:12,:);
            cf = dir.*obj.CF(1:6,:) + (1 - dir).*obj.CF(7:12,:);
            ma = dir.*obj.MA(1:6,:) + (1 - dir).*obj.MA(7:12,:);
            
            % Calculate Hydrostatic forces
            Tg=  [-sin(obj.n(5))*obj.W;
                   cos(obj.n(5))*sin(obj.n(4))*obj.W;
                   cos(obj.n(5))*cos(obj.n(4))*obj.W;
                   0;
                   0;
                   0]; %--> Weight force/moment vector              
            Tb=  [sin(obj.n(5))*obj.B;
                  -cos(obj.n(5))*sin(obj.n(4))*obj.B;
                  -cos(obj.n(5))*cos(obj.n(4))*obj.B;
                  -obj.COB(2)*cos(obj.n(5))*cos(obj.n(4))*obj.B + obj.COB(3)*cos(obj.n(5))*sin(obj.n(4))*obj.B;
                  obj.COB(3)*sin(obj.n(5))*obj.B + obj.COB(1)*cos(obj.n(5))*cos(obj.n(4))*obj.B;
                  -obj.COB(1)*cos(obj.n(5))*sin(obj.n(4))*obj.B - obj.COB(2)*sin(obj.n(5))*obj.B]; %Buoyancy force/moment vector
            T_Static = Tg + Tb; %--> Total Hydrostatic Force/Moment
            
            % Calculate Hydrodynamic forces
            Td = -cd * (obj.v.'.*abs(obj.v.')); %--> Drag force/moment vector     
            Tl = cf * (obj.v.'.*obj.v.'); %--> Lift force/moment vector
            T_Dynamic = Td + Tl; %--> Total Hydrodynamic Force/Moment
            
            % Write Thrust forces
            obj = ESC(obj,PWM);
            Tt = obj.proj*getTh(obj).'; %--> Resultant thrust force/moment
            
            % Newton's 2nd law
            T = T_Static + T_Dynamic + Tt; %--> Total acting forces/moments
            Crb = [zeros(3,3) (-obj.m*skew(obj.v(1:3).'));  
                   (-obj.m*skew(obj.v(1:3).')) (-skew(obj.I0*obj.v(4:6).'))]; %--> Coriolis
            Ca = [zeros(3,3), (-skew(ma(1:3,1:3)*obj.v(1:3).'+ma(1:3,4:6)*obj.v(4:6).'));  
                 (-skew(ma(1:3,1:3)*obj.v(1:3).'+ma(1:3,4:6)*obj.v(4:6).')),...
                 (-skew(ma(4:6,1:3)*obj.v(1:3).'+ma(4:6,4:6)*obj.v(4:6).'))]; %--> Added Mass Coriolos
            obj.vdot = ((obj.M + ma) \ (T - (Ca+Crb)*obj.v.')).'; %--> AUV acceleration
            obj.ndot = (Jacobian(obj)*obj.v.').'; %--> Derivative of global states
        end
        function [obj,dT_new] = Solve(obj,PWM,dT)
            % Ranga-kutta parameters 
            dn_max = 0.001; %--> max relative change tolerance
            dn_min = 0.00001; %--> min relative change tolerance
            
            % Save starting states 
            n_start = obj.n; %--> Starting global pose
            v_start = obj.v; %--> Starting local velocity
            
            % Normal step 
            obj = Newton(obj,PWM); %--> Get acceleration (K1)
            ndot_1 = obj.ndot; vdot_1 = obj.vdot; %--> Save starting acceleration (K1)
            avgdn = ndot_1; avgdv = vdot_1; %--> Averaging variable to get average of accelerations (K1,2K2,2K3,k4)
            obj.n = n_start + (dT/2)*obj.ndot; obj.v = v_start + (dT/2)*obj.vdot; %--> 1/2 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K2)
            avgdn = avgdn + 2*obj.ndot; avgdv = avgdv + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + (dT/2)*obj.ndot; obj.v = v_start + (dT/2)*obj.vdot; %--> 1/2 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K3)
            avgdn = avgdn + 2*obj.ndot; avgdv = avgdv + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + dT*obj.ndot; obj.v = v_start + dT*obj.vdot; %--> 1 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K4)
            avgdn = avgdn + obj.ndot; avgdv = avgdv + obj.vdot; %--> Update averaging variable
            step_n = n_start + (dT/6)*avgdn; %--> Compute full step global pose update
            step_v = v_start + (dT/6)*avgdv; %--> Compute full step local velocity update
            
            % Half step 
            obj.ndot = ndot_1; obj.vdot = vdot_1; %--> reset acceleration value to K1
            obj.n = n_start + (dT/4)*obj.ndot; obj.v = v_start + (dT/4)*obj.vdot; %--> 1/4 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K2)
            avgdn = ndot_1 + 2*obj.ndot; avgdv = vdot_1 + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + (dT/4)*obj.ndot; obj.v = v_start + (dT/4)*obj.vdot; %--> 1/4 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K3)
            avgdn = avgdn + 2*obj.ndot; avgdv = avgdv + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + (dT/2)*obj.ndot; obj.v = v_start + (dT/2)*obj.vdot; %--> 1/2 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K3)
            avgdn = avgdn + obj.ndot; avgdv = avgdv + obj.vdot; %--> Update averaging variable
            half_step_n = n_start + (dT/12)*avgdn; %--> Compute half step global pose update
            half_step_v = v_start + (dT/12)*avgdv; %--> Compute half step local velocity update
            
            % Double step 
            obj.ndot = ndot_1; obj.vdot = vdot_1; %--> reset acceleration value to K1
            obj.n = n_start + dT*obj.ndot; obj.v = v_start + dT*obj.vdot; %--> 1 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K2)
            avgdn = ndot_1 + 2*obj.ndot; avgdv = vdot_1 + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + dT*obj.ndot; obj.v = v_start + dT*obj.vdot; %--> 1 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K3)
            avgdn = avgdn + 2*obj.ndot; avgdv = avgdv + 2*obj.vdot; %--> Update averaging variable
            obj.n = n_start + (dT*2)*obj.ndot; obj.v = v_start + (dT*2)*obj.vdot; %--> 2 sampling time development
            obj = Newton(obj,PWM); %--> Get acceleration (K4)
            avgdn = avgdn + obj.ndot; avgdv = avgdv + obj.vdot; %--> Update averaging variable
            dbl_step_n = n_start + (dT/3)*avgdn; %--> Compute double step global pose update
            dbl_step_v = v_start + (dT/3)*avgdv; %--> Compute double step local velocity update
            
            % Correct step size 
            if(dT<obj.dT_min) %--> Chech sampling time doesn't fall below a minimum value
                dT_new = obj.dT_min;
                obj.n = step_n; obj.v = step_v;
                obj.time = obj.time + dT;
            else
                if(min(abs(step_n-half_step_n)) > dn_max) %--> Use half step if error exceeds threshold (dn_max)
                    dT_new = dT/2;
                    obj.n = half_step_n; obj.v = half_step_v;
                    obj.time = obj.time + dT_new;
                else
                    if(max(abs(step_n-dbl_step_n)) < dn_min) %--> Use dbl step if error is less than threshold (dn_min)
                        dT_new = dT*2;
                        obj.n = dbl_step_n; obj.v = dbl_step_v;
                        obj.time = obj.time + dT_new;
                    else
                        obj.n = step_n; obj.v = step_v; %--> Stick to normal step if you lie between both thresholds
                        dT_new = dT;
                        obj.time = obj.time + dT_new;
                    end
                end
            end
            
            % Update AUV clock 
            obj.iter = obj.iter + 1; %--> Update number of iterations              
        end
    end
end

% Ordinary functions
function S = skew(a)
    % Return a skew symmetric matrix made up of elements of vector a
    S = [0 -a(3,1) a(2,1); 
        a(3,1) 0 -a(1,1); 
        -a(2,1) a(1,1) 0];
end