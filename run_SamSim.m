          % Author: Samer A. Mohamed (https://github.com/SamMans) %
%-------------------------------------------------------------------------%
%----------------------------  AUV Simulation  ---------------------------%
%-------------------------------------------------------------------------%
% This script is used to run the AUV simulation environment
% The parameters of the state space model depend on choices from SamSim.mlx

% Load AUV model parameters
load(pwd + "/Libraries/AUVs/" + model + "/physical.mat");

% Initial conditions
dT = 10^-1; %--> Initial sampling time of the Runga-Kutta algorithm
n = zeros(1, 6); %--> Initial global pose
v = zeros(1, 6); %--> Initial local velocity
v_dot = zeros(1, 6); %--> Initial acceleration
Ti = 1; %--> Initial mission planner target index
Target = zeros(size(MP_margin, 1), 4); %--> Initial target

% Create objects
Vehicle = AUV(model, Landmarks, n, v, v_dot); %--> AUV object
if(nav_pkg == "Bolt_1")
    Vehicle_nav = Nav_Stack_2020(model, Th_lim, getn(Vehicle), Landmarks); %--> Navigation stack object
    Kin = zeros((Sim_Time / getdT_min(Vehicle)) + 1, 3); %--> Kinematic controls history
elseif(nav_pkg == "Bolt_2")
    Vehicle_nav = Nav_Stack_2021(model, Th_lim, getn(Vehicle), Landmarks); %--> Navigation stack object
end

% Histories (for records, plots and logs)
NN = zeros((Sim_Time / getdT_min(Vehicle)) + 1, 7); %--> Global pose history
NN(1, :) = [0, n]; %--> Log initial pose into history vector
VV = zeros((Sim_Time / getdT_min(Vehicle)) + 1, 7); %--> Velocity history
VV(1, :) = [0, v]; %--> Log initial velocity into history vector
TT = zeros((Sim_Time / getdT_min(Vehicle)), 9); %--> Thrust history
SS = zeros((Sim_Time / getdT_min(Vehicle)), 5); %--> Sensor history
Bel = zeros((Sim_Time / getdT_min(Vehicle)) + 1, 9 + ...
    size(Landmarks, 1) * size(Landmarks, 2)); %--> SLAM belief history 
Bel(1, :) = [0, v(1 : 3), v(6), n(1 : 3), n(6), ...
    zeros(1, size(Landmarks, 1) * size(Landmarks, 2))]; %--> Log initial belief into history vector
Cov = zeros((Sim_Time / getdT_min(Vehicle)) + 1, 9 + ...
    size(Landmarks, 1) * size(Landmarks, 2)); %--> SLAM covariance history
Cov(1, :) = [0, zeros(1, 8), 10000 * ones(1, size(Landmarks, 1) * ...
    size(Landmarks, 2))]; %--> Log initial covariance into history vector

% Run model
while(gett(Vehicle) < Sim_Time)
    % Get sensor readings at time "t_now"
    [Zread, RPYread, LMread] = Refresh(Vehicle.Sens, getn(Vehicle), gett(Vehicle));
    t_now = gett(Vehicle);
    
    % Generate AUV pwm signals "from control node"
    Vehicle_nav = Control(Vehicle_nav, gett(Vehicle), Target(Ti, :), Zread, RPYread, LMread);
                  
    % Update AUV states
    [Vehicle, dT] = Solve(Vehicle, getpwm(Vehicle_nav), dT);
                           
    % Update history of states 
    NN(geti(Vehicle) + 1, :) = [gett(Vehicle), getn(Vehicle)]; 
    VV(geti(Vehicle) + 1, :) = [gett(Vehicle), getv(Vehicle)];
    TT(geti(Vehicle), :) = [t_now, getTh(Vehicle)];
    SS(geti(Vehicle), :) = [t_now, Zread, RPYread];
    Bel(geti(Vehicle) + 1, :) = [gett(Vehicle), getbel(Vehicle_nav)];
    Cov(geti(Vehicle) + 1, :) = [gett(Vehicle), getcov(Vehicle_nav)];
    if(nav_pkg == "Bolt_1")
        Kin(geti(Vehicle) + 1, :) = [gett(Vehicle), getkincon(Vehicle_nav)];
    end

    % Mission planner activation
    if(Bel(geti(Vehicle) + 1, 6) <= MP_margin(Ti, 1) + Target(Ti, 1) && ...
            Bel(geti(Vehicle) + 1, 6) >= -MP_margin(Ti, 1) + Target(Ti, 1) && ...
            Bel(geti(Vehicle) + 1, 7) <= MP_margin(Ti, 2) + Target(Ti, 2) && ...
            Bel(geti(Vehicle) + 1, 7) >= -MP_margin(Ti, 2) + Target(Ti, 2) && ...
            Bel(geti(Vehicle) + 1, 8) <= MP_margin(Ti, 3) + Target(Ti, 3) && ...
            Bel(geti(Vehicle) + 1, 8) >= -MP_margin(Ti, 3) + Target(Ti, 3) && ...
            Bel(geti(Vehicle) + 1, 9) <= MP_margin(Ti, 4) + Target(Ti, 4) && ...
            Bel(geti(Vehicle) + 1, 9) >= -MP_margin(Ti, 4) + Target(Ti, 4))
        if(Ti < size(MP_margin, 1))
            if(~(Target(Ti, 1) == 0 && Target(Ti, 2) == 0) || ...
                    (Target(Ti, 1) == 0 && Target(Ti, 2) == 0 && Ti ~= 1))
                Ti = Ti + 1; %--> Move out to next goal, if old goal is achieved
            end
        else
            disp("Mission Success!!"); %--> Mission termination
            break;
        end
    end
    XYZ = reshape(Bel(geti(Vehicle) + 1, 10 : size(Bel, 2)), [], 3);
    Target = MP(XYZ(1, :), XYZ(2, :), XYZ(3, :)); %--> Upcoming target
end

% Post processing of results
% Cropping extra history vector length
NN = NN(1 : geti(Vehicle) + 1, :); VV = VV(1 : geti(Vehicle) + 1, :); 
TT = TT(1 : geti(Vehicle), :); SS = SS(1 : geti(Vehicle), :);
Bel = Bel(1 : geti(Vehicle) + 1, :); Cov = Cov(1 : geti(Vehicle) + 1, :);
if(nav_pkg == "Bolt_1")
    Kin = Kin(1 : geti(Vehicle) + 1, :);
end

% Additional options
% Excel logs
if(Logs == true)
    % Ground truth pose record
    log_t = string(clock);
    writematrix(["Time (sec.)", "X (m)", "Y (m)", "Z (m)", "roll (rad)", "pitch (rad)", "yaw (rad)"], ...
        pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Ground_Truth_Pose');
    writematrix(NN, pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Ground_Truth_Pose');

    % Ground truth velocity record
    writematrix(["Time (sec.)", "U (m/s)", "V (m/s)", "W (m/s)", "P (rad/s)", ...
        "Q (rad/s)", "R (rad/s)"], ...
        pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Ground_Truth_Vel');
    writematrix(VV, pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Ground_Truth_Vel');
    
    % Thrust record
    writematrix(["Time (sec.)", "T1 (N)", "T2 (N)", "T3 (N)", "T4 (N)", ...
        "T5 (N)", "T6 (N)", "T7 (N)", "T8 (N)"], ...
        pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Thrust');
    writematrix(TT, pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Thrust');

    % Sensors record
    writematrix(["Time (sec.)", "Depth (m)", "roll (rad)", "pitch (rad)", "yaw (rad)"], ...
        pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Sensors');
    writematrix(SS, pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Sensors');

    % Belief record
    writematrix(["Time (sec.)", "U (m/sec.)", "V (m/s)", "W (m/s)", "R (rad/s)", "X (m)", ...
        "Y (m)", "Z (m)", "Yaw (rad)"], pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Belief');
    writematrix(Bel(:, 1 : 9), pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Belief');

    % Covariance record
    writematrix(["Time (sec.)", "U", "V", "W", "R", "X", ...
        "Y", "Z", "Yaw"], pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'Sheet', 'Covariance');
    writematrix(Cov(:, 1 : 9), pwd + "/Log Files/" + log_t(4) + "_" + log_t(5) + "_" + ...
        log_t(6) + "_" + log_t(3) + "_" + log_t(2) + "_" + log_t(1) + ".xlsx", ...
        'WriteMode', 'append', 'Sheet', 'Covariance');
end

% Graphs
Sz = 330; %--> Plot size
if(Plots == true)
    % SLAM plots
    % Pose (Ground truth and estimated) vs. time
    Pose_plot = ["X (m)", "Y (m)", "Z (m)", "Roll (rad)", "Pitch (rad)", "Yaw (rad)"];
    Est = [Bel(:, 1), Bel(:, 6), Bel(:, 7), Bel(:, 8), zeros(size(Bel, 1), 1), ...
        zeros(size(Bel, 1), 1), Bel(:, 9)];
    for i = 1 : 6
        figure
        plot(NN(:, 1), NN(:, i + 1), 'r', 'LineWidth', 2)
        hold on 
        plot(Est(:, 1), Est(:, i + 1), '--b', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel(Pose_plot(1, i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        legend('Ground truth','Estimate', 'Location', 'Best')
        set(gcf,'color','w');
    end

    % Velocity (Ground truth and estimated) vs. time
    Vel_plot = ["Surge (m/s)", "Starboard (m/s)", "Heave (m/s)", ...
    "About surge (rad/s)", "About starboard (rad/s)", "About heave (rad/s)"];
    Est = [Bel(:, 1), Bel(:, 2), Bel(:, 3), Bel(:, 4), zeros(size(Bel, 1), 1), ...
        zeros(size(Bel, 1), 1), Bel(:, 5)];
    for i = 1 : 6
        figure
        plot(VV(:, 1), VV(:, i + 1), 'r', 'LineWidth', 2)
        hold on 
        plot(Est(:, 1), Est(:, i + 1), '--g', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel(Vel_plot(1, i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        legend('Ground truth','Estimate', 'Location', 'Best')
        set(gcf,'color','w');
    end

    % Mission location estimation vs. time
    for i = 10 : 3 : size(Bel, 2)
        figure
        plot(Bel(:, 1), Bel(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 10)/3 + 1) + " X (m)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
    for i = 11 : 3 : size(Bel, 2)
        figure
        plot(Bel(:, 1), Bel(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 11)/3 + 1) + " Y (m)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
    for i = 12 : 3 : size(Bel, 2)
        figure
        plot(Bel(:, 1), Bel(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 12)/3 + 1) + " Z (m)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end

    % Covariance vs. time
    Cov_plot = ["U (m2/s2)", "V (m2/s2)", "W (m2/s2)", "R (rad2/s2)", ...
        "X (m2)", "Y (m2)", "Z (m2)", "Yaw (rad2)"];
    for i = 1 : 8
        figure
        plot(Cov(:, 1), Cov(:, i + 1), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel(Cov_plot(1, i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
    for i = 10 : 3 : size(Bel, 2)
        figure
        plot(Cov(:, 1), Cov(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 10)/3 + 1) + " X (m2)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
    for i = 11 : 3 : size(Bel, 2)
        figure
        plot(Cov(:, 1), Cov(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 11)/3 + 1) + " Y (m2)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
    for i = 12 : 3 : size(Bel, 2)
        figure
        plot(Cov(:, 1), Cov(:, i), 'k', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 10), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("Landmark " + string((i - 12)/3 + 1) + " Z (m2)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end

    % Controls vs. time
    for i = 1 : 8
        figure
        plot(TT(:, 1), TT(:, i + 1), 'r', 'LineWidth', 2)
        annotation('textbox', [.0 .0 .0 1.0], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
        xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        ylabel("T" + string(i) + " (N)", 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
        set (gca, 'fontweight', 'bold', 'FontSize', 18) 
        set(gcf,'units','points','position',[.0, .0, Sz, Sz])
        set(gcf,'color','w');
    end
end

% Run simulink
if(Video == true)
    % Video duration
    dur_vid = gett(Vehicle);

    % State feed to visualization 
    NNx = [NN(:, 1), NN(:, 2)]; NNy = [NN(:, 1), NN(:, 3)];
    NNz = [NN(:, 1), NN(:, 4)]; NNphi = [NN(:, 1), NN(:, 5)];
    NNtht = [NN(:, 1), NN(:, 6)]; NNeps = [NN(:, 1), NN(:, 7)];

    % Run simulation
    sim(pwd + "/Pools/" + pool + "/test_field.slx");
end
