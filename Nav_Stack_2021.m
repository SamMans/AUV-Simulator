         % Author: Samer A. Mohamed (https://github.com/SamMans) %
classdef Nav_Stack_2021
    properties(SetAccess = private)
    % Vehicle parameters
        B_table %--> Mapping vector (thruster <--> pwm)
        th_lim %--> Thrust limit per thruster "software limit"
        
    % Controller parameters
        dT_con %--> Controller sampling time
        Ks %--> Stabilizer gain matrix
        RPYread_k1 %--> Past IMU readings
        pwm_th %--> Output pwm to ESCs plus actual thrust
        solver %--> Non-linear MPC solver
        N_hor %--> MPC prediction horizon
        Mx %--> Inertia in surge direction
        My %--> Inertia in starboard direction
        Mth %--> Inertia about yaw
        Dx %--> Drag in surge direction
        Dy %--> Drag in starboard direction
        Dth %--> Drag about yaw
        tht %--> Thruster mounting angle
        rx %--> Vertical thruster arm in surge direction
        ry %--> Vertical thruster arm in starboard direction
        rh %--> Horizontal thrustr arm length
    
    % SLAM parameters 
        u %--> Belief mean "current"
        Sigma %--> Belief covariance
        Dz %--> Drag in heave direction
        Mz %--> Inertia in heave direction
        m %--> AUV mass
        V %--> AUV volume
        T_com %--> Transform from camera frame to com frame
        TRANS_p %--> Pressure sensor shift from center of mass
        SensVar %--> Variance of sensor readings
    
    % Timing parameters
        Tr %--> Last recorded system time
    end
    methods
        function obj = Nav_Stack_2021(model, Th_lim, n, lm) %---> Constructor
            % Load vehicle parameters
            load(pwd + "\Libraries\AUVs\" + model + "\physical.mat");
            obj.B_table = readmatrix(pwd + "\Libraries\Thrusters\" + ...
                Th_mod(1) + ".xls", 'Sheet', Th_mod(2));
            obj.th_lim = Th_lim;
            obj.pwm_th = [1500 * ones(1, 8); zeros(1, 8)];
            
            % Timing parameters
            obj.Tr = 0;
            
            % EKF & sensor reading parameters initializations
            obj.u = [0, 0, 0, 0, n(1), n(2), n(3), n(6), zeros(1, 3 * size(lm, 1))];
            obj.Sigma = 10000 * eye(size(obj.u, 2)); 
            for i = 1:8
                obj.Sigma(i, i) = 0;
            end
            obj.RPYread_k1 = [n(4), n(5), n(6)];
            obj.Dz = [CD_AUV(3,3), CD_AUV(9,3)]; 
            obj.Mz = [m_AUV + MA_AUV(3,3), m_AUV + MA_AUV(9,3)];
            obj.m = m_AUV; obj.V = V_AUV;
            obj.T_com = [1 0 0 Sshift(1);         
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
                     0 0 0 1];
            load(pwd + "\Libraries\Sensors\Depth sensors\" + Sens_mod(1) + ".mat");
            load(pwd + "\Libraries\Sensors\IMUs\" + Sens_mod(2) + ".mat");
            load(pwd + "\Libraries\Sensors\Stereo cameras\" + Sens_mod(3) + ".mat");
            obj.SensVar = [(Pprec / 300)^2, (imu_dev)^2, (prec_s / 300)^2];
            obj.TRANS_p = Pshift(3);
            
            % Stabilizer parameters setup
            obj.dT_con = 0.1;
            Kp_z = 1000; Kd_z = 2*sqrt(Kp_z * (MA_AUV(3,3) + m_AUV));
            Kp_r = 1000; Kd_r = 2*sqrt(Kp_r * (MA_AUV(4,4) + I0_AUV(1,1)));
            Kp_p = 1000; Kd_p = 2*sqrt(Kp_p * (MA_AUV(5,5) + I0_AUV(2,2))); 
            obj.Ks = [-1, (Kp_z + (2*Kd_z) / obj.dT_con), (Kp_z - (2*Kd_z) / obj.dT_con);
                -1, (Kp_r + (2*Kd_r) / obj.dT_con), (Kp_r - (2*Kd_r) / obj.dT_con);
                -1, (Kp_p + (2*Kd_p) / obj.dT_con), (Kp_p - (2*Kd_p) / obj.dT_con)];    
            
            % Trajectory tracking parameter setup
            obj.Mx = m_AUV + MA_AUV(1,1); obj.My = m_AUV + MA_AUV(2,2);
            obj.Mth = I0_AUV(3,3) + MA_AUV(6,6);
            obj.Dx = CD_AUV(1,1); obj.Dy = CD_AUV(2,2); obj.Dth = CD_AUV(6,6); 
            obj.tht = Tht_AUV;
            obj.rx = Rx_AUV; obj.ry = Ry_AUV; obj.rh = R_Hor_AUV;
            N = 5; %--> MPC prediction horizon
            obj.N_hor = N;
            Q = 5000000 * eye(3,3); %--> Weighing matrices (states)
            R = eye(4,4); %--> Weighing matrices (controls)
            
            % Trajectory tracking solver setup
            addpath(pwd + "\Dependencies\Casadi")
            import casadi.*
            x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
            vx = SX.sym('vx'); vy = SX.sym('vy'); w = SX.sym('w');
            states = [x; y; theta; vx; vy; w]; n_states = length(states); %--> Symbolic states
            T1 = SX.sym('T1'); T2 = SX.sym('T2'); T3 = SX.sym('T3'); T4 = SX.sym('T4');
            controls = [T1; T2; T3; T4]; n_controls = length(controls); %--> Symbolic controls
            rhs = [vx * cos(theta) - vy * sin(theta);
              vx * sin(theta) + vy * cos(theta);
              w;
              (-obj.Dx * vx * abs(vx) + cos(obj.tht) * (T1 + T2 + T3 + T4)) / obj.Mx;
              (-obj.Dy * vy * abs(vy) + cos(obj.tht) * (-T1 + T2 - T3 + T4)) / obj.My;
              (-obj.Dth * w * abs(w) + obj.rh * (-T1 + T2 + T3 - T4)) / obj.Mth]; %--> State space model
            f = Function('f', {states, controls}, {rhs}); %--> Nonlinear mapping function f(x,u)
            U = SX.sym('U', n_controls, N); %--> Decision variables (controls)
            P = SX.sym('P', n_states + 3); %--> Initial state and reference state
            X = SX.sym('X', n_states, (N + 1)); %--> States over the optimization problem
            PHI = 0; %--> Objective function
            g = [];  %--> Constraints vector
            st  = X(:, 1); %--> Initial state
            g = [g; st - P(1:n_states)]; %--> Initial condition constraints
            for k = 1:N
                st = X(:, k);  con = U(:, k);
                PHI = PHI + (st(1:3) - P(7:9))' * Q * (st(1:3) - P(7:9)) + ...
                    con' * R * con; %--> Calculate obj
                st_next = X(:, k+1);
                f_value = f(st, con);
                st_next_euler = st + (obj.dT_con * f_value);
                g = [g; st_next - st_next_euler]; % Compute constraints
            end
            OPT_variables = [reshape(X, n_states * (N+1), 1); reshape(U, n_controls * N, 1)];
            nlp_prob = struct('f', PHI, 'x', OPT_variables, 'g', g, 'p', P);
            opts = struct;
            opts.ipopt.max_iter = 2000;
            opts.ipopt.print_level = 0;
            opts.print_time = 0;
            opts.ipopt.acceptable_tol = 1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
            obj.solver = nlpsol('solver', 'ipopt', nlp_prob, opts);
        end
        function obj = T2P(obj, Thrus) %--> Thrust to pwm converter + thrust correction 
            % Limit thrust
            Thrus(Thrus > obj.th_lim) = obj.th_lim; 
            Thrus(Thrus < -obj.th_lim) = -obj.th_lim;
            
            % Loop to fill the pwm vector
            for i = 1 : size(obj.pwm_th, 2)
                [~,j] = min(abs(Thrus(i) - obj.B_table(:,6) * 9.80665));
                obj.pwm_th(1, i) = obj.B_table(j, 1);
                obj.pwm_th(2, i) = obj.B_table(j, 6)* 9.80665;
            end
        end
        function PWM = getpwm(obj) %--> Get control pwm from navigation stack
            PWM = obj.pwm_th(1, :);
        end
        function bel = getbel(obj) %--> Get current belief
            bel = obj.u;
        end
        function cov = getcov(obj) %--> Get current covariance
            cov = diag(obj.Sigma).'; 
        end
        function obj = Control(obj, Rt, Setp, zread, rpyread, lmread)
            if(round(Rt - obj.Tr, 1) >= obj.dT_con)
                % SLAM call
                Th_k1 = obj.pwm_th(2, :); %--> Last thrust
                u_k1 = obj.u; %--> Last belief
                for i = 1 : size(lmread, 1)
                    transf = obj.T_com * [lmread(i, 2:4), 1].'; %--> Apply transformation
                    lmread(i, 2:4) = transf(1:3, 1).'; %--> Update camera readings to be referenced to COM
                end
                obj = EKF(obj, Th_k1, zread, rpyread, lmread); %--> Belief update
                
                % Compute errors
                ez = Setp(3) - obj.u(7); %--> Current depth error
                ez_k1 = Setp(3) - u_k1(7); %--> Past depth error
                er = -rpyread(1); %--> Current roll error
                er_k1 = -obj.RPYread_k1(1); %--> Past roll error
                ep = -rpyread(2); %--> Current pitch error
                ep_k1 = -obj.RPYread_k1(2); %--> Past pitch error
                
                % Apply PID difference equation
                Ez = [Th_k1(5) + Th_k1(6) + Th_k1(7) + ...
                       Th_k1(8); ez; ez_k1];
                Er = [(Th_k1(5) - Th_k1(6) - Th_k1(7) + ...
                       Th_k1(8)) * obj.ry; er; er_k1];
                Ep = [(-Th_k1(5) - Th_k1(6) + Th_k1(7) + ...
                       Th_k1(8)) * obj.rx; ep; ep_k1];
                U = [obj.Ks(1,:) * Ez, obj.Ks(2,:) * Er, obj.Ks(3,:) * Ep];

                % Assign thrust values
                Th_ver = (U(1) / 4) * ones(1,4) + ...
                    (U(2) / (4 * obj.ry)) * [1, -1, -1, 1] + ...
                    (U(3) / (4 * obj.rx)) * [-1, -1, 1, 1];
                
                % Nonlinear MPC kinematic control
                n_states = 6; n_controls = 4;
                args = struct;
                args.lbg(1:n_states*(obj.N_hor+1)) = 0; %--> Equality constraints
                args.ubg(1:n_states*(obj.N_hor+1)) = 0; %--> Equality constraints
                args.lbx(1:n_states*(obj.N_hor+1),1) = -inf; %--> States lower bound
                args.ubx(1:n_states*(obj.N_hor+1),1) = inf; %--> States upper bound
                args.lbx(n_states*(obj.N_hor+1)+1:n_states*(obj.N_hor+1)+ ...
                    n_controls*obj.N_hor,1) = -obj.th_lim; %--> Thrust lower bound
                args.ubx(n_states*(obj.N_hor+1)+1:n_states*(obj.N_hor+1)+ ...
                    n_controls*obj.N_hor,1) = obj.th_lim; %--> Thrust upper bound
                x0 = [obj.u(5) ; obj.u(6) ; obj.u(8); obj.u(1); obj.u(2); obj.u(4)]; %--> Initial condition
                xs = [Setp(1); Setp(2); Setp(4)]; %--> Reference posture
                u0 = zeros(obj.N_hor,n_controls); %--> Four control inputs
                X0 = repmat(x0,1,obj.N_hor+1)'; %--> Initialization of the decision variables
                args.p = [x0; xs]; %--> Set the values of the parameters vector
                args.x0  = [reshape(X0',n_states*(obj.N_hor+1),1);reshape(u0',n_controls*obj.N_hor,1)];
                sol = obj.solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
                    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
                Th_hor = reshape(full(sol.x(n_states*(obj.N_hor+1)+1:end))',n_controls,obj.N_hor)'; %--> Get controls only from the solution
                Th_hor = Th_hor(1,:); %--> Horizontal thruster command
                
                % Assign full eight thrust values
                obj = T2P(obj, [Th_hor, Th_ver]);
                
                % Update private variables
                obj.RPYread_k1 = rpyread;
                obj.Tr = Rt;
            end
        end
        function obj = EKF(obj, Th, depthr, rpyr, lmr)
            % Prediction step
            % 1: Update covariance
            F = zeros(8, size(obj.u, 2)); %--> F simplification matrix 
            R = 10^-4 * eye(8); %--> Noise matrix
            for i = 1:8
                F(i, i) = 1;
            end
            jac = [(-2 * obj.dT_con * obj.Dx * abs(obj.u(1))) / obj.Mx, 0, 0, 0, 0, 0, 0, 0;
                0, (-2 * obj.dT_con * obj.Dy * abs(obj.u(2))) / obj.My, 0, 0, 0, 0, 0, 0;
                0, 0, (-2 * obj.dT_con * obj.Dz(round((-sign(obj.u(3))+2)/2)) * abs(obj.u(3))) / obj.Mz(round((-sign(obj.u(3))+2)/2)), 0, 0, 0, 0, 0;
                0, 0, 0, (-2 * obj.dT_con * obj.Dth * abs(obj.u(4))) / obj.Mth, 0, 0, 0, 0;
                cos(obj.u(8)) * obj.dT_con, -sin(obj.u(8)) * obj.dT_con, 0, 0, 0, 0, 0, -(obj.u(1) * sin(obj.u(8)) + obj.u(2) * cos(obj.u(8))) * obj.dT_con;
                sin(obj.u(8)) * obj.dT_con, cos(obj.u(8)) * obj.dT_con, 0, 0, 0, 0, 0, (obj.u(1) * cos(obj.u(8)) - obj.u(2) * sin(obj.u(8))) * obj.dT_con;
                0, 0, obj.dT_con, 0, 0, 0, 0, 0;
                0, 0, 0, obj.dT_con, 0, 0, 0, 0]; %--> Dynamic model jacobian
            G = eye(size(obj.u, 2)) + F.' * jac * F; %--> Full jacobian
            obj.Sigma = G * obj.Sigma * G.' + F.' * R * F; %--> Belief covariance update
            
            % 2: Motion model
            W = obj.m * 9.81; B = 1000 * obj.V * 9.81; %--> Static forces
            Taw = [sum(Th(1:4)) * cos(obj.tht), sum(Th(1:4).*[-1, 1, -1, 1]) * cos(obj.tht), ...
                sum(Th(5:8)), sum(Th(5:8).*[1, -1, -1, 1]) * obj.ry, sum(Th(5:8).*[-1, -1, 1, 1]) * obj.rx, ...
                sum(Th(1:4).*[-1, 1, 1, -1]) * obj.rh]; %--> External forces
            obj.u = obj.u + [(Taw(1) - obj.Dx * obj.u(1) * abs(obj.u(1))) * (obj.dT_con / obj.Mx), ...
                (Taw(2) - obj.Dy * obj.u(2) * abs(obj.u(2))) * (obj.dT_con / obj.My), ...
                (Taw(3) - obj.Dz(round((-sign(obj.u(3))+2)/2)) * obj.u(3) * abs(obj.u(3)) + W - B) ...
                * (obj.dT_con / obj.Mz(round((-sign(obj.u(3))+2)/2))), ...
                (Taw(6) - obj.Dth * obj.u(4) * abs(obj.u(4))) * (obj.dT_con / obj.Mth), ...
                (obj.u(1) * cos(obj.u(8)) - obj.u(2) * sin(obj.u(8))) * obj.dT_con, ...
                (obj.u(1) * sin(obj.u(8)) + obj.u(2) * cos(obj.u(8))) * obj.dT_con, ...
                obj.u(3) * obj.dT_con, ...
                obj.u(4) * obj.dT_con]* F; %--> Belief mean update
            
            % Update step
            % Step 1: key variables: z, h, H, Q
            z = zeros(2 + size(lmr,1) * 3, 1); z(1,1) = depthr; z(2,1) = rpyr(3); %--> Measurements vector
            h = zeros(2 + size(lmr,1) * 3, 1); h(1,1) = obj.u(7) + obj.TRANS_p; h(2,1) = obj.u(8); %--> Measurements prediction vector
            H = zeros(2 + size(lmr,1) * 3, size(obj.u, 2)); %--> Measurement jacobian
            H(1, 7) = 1; H(2, 8) = 1;  
            Q = zeros(2 + size(lmr,1) * 3, 2 + size(lmr,1) * 3); %--> Sensor noise matrix
            Q(1, 1) = obj.SensVar(1); Q(2, 2) = obj.SensVar(2);
            for i = 3 : size(Q, 1)
                Q(i, i) = obj.SensVar(3);
            end
            for i = 1 : size(lmr,1)
                % Update sensor readings using camera input
                z(3 + (i-1) * 3, 1) = lmr(i, 2); z(4 + (i-1) * 3, 1) = lmr(i, 3);
                z(5 + (i-1) * 3, 1) = lmr(i, 4);
                
                % Update measurements vector according to observations
                mx = obj.u(9 + 3 * lmr(i,1)); my = obj.u(10 + 3 * lmr(i,1)); 
                mz = obj.u(11 + 3 * lmr(i,1));
                h(3 + (i-1) * 3, 1) = mx * cos(obj.u(8)) + my * sin(obj.u(8)) - obj.u(5) * cos(obj.u(8)) ...
                    - obj.u(6) * sin(obj.u(8));
                h(4 + (i-1) * 3, 1) = -mx * sin(obj.u(8)) + my * cos(obj.u(8)) + obj.u(5) * sin(obj.u(8))...
                    - obj.u(6) * cos(obj.u(8));
                h(5 + (i-1) * 3, 1) = mz - obj.u(7);
               
                % Update measurements jacobian according to observations
                H(3 + (i-1) * 3, 5) = -cos(obj.u(8)); H(3 + (i-1) * 3, 6) = -sin(obj.u(8));
                H(3 + (i-1) * 3, 8) = -mx*sin(obj.u(8)) + my*cos(obj.u(8)) + obj.u(5)*sin(obj.u(8)) - obj.u(6)*cos(obj.u(8));
                H(3 + (i-1) * 3, 9 + 3 * lmr(i,1)) = cos(obj.u(8)); H(3 + (i-1) * 3, 10 + 3 * lmr(i,1)) = sin(obj.u(8));
                H(4 + (i-1) * 3, 5) = sin(obj.u(8)); H(4 + (i-1) * 3, 6) = -cos(obj.u(8)); 
                H(4 + (i-1) * 3, 8) = -mx*cos(obj.u(8)) - my*sin(obj.u(8)) + obj.u(5)*cos(obj.u(8)) + obj.u(6)*sin(obj.u(8));
                H(4 + (i-1) * 3, 9 + 3*lmr(i,1)) = -sin(obj.u(8)); H(4 + (i-1) * 3, 10 + 3*lmr(i,1)) = cos(obj.u(8));
                H(5 + (i-1) * 3, 7) = -1; H(5 + (i-1) * 3, 11 + 3*lmr(i,1)) = 1;
            end
            
            % Step 2: Kalman equations & final belief update
            K = obj.Sigma * H.' * inv(H * obj.Sigma * H.' + Q); %--> Kalman gain
            obj.u = (obj.u.' + K * (z - h)).';
            obj.Sigma = (eye(size(obj.u, 2)) - K * H) * obj.Sigma;
        end
    end
end