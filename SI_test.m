          % Author: Samer A. Mohamed (https://github.com/SamMans) %
%-------------------------------------------------------------------------%
%---------------------------  AUV SI Simulation  -------------------------%
%-------------------------------------------------------------------------%
% This script is used to run the AUV simulation environment
% This model only tests SI success by comparing against actual Excel logs
% The test data time stamps have to be equally spaced at 0.1 sec. intervals

% User choices
model = "Bolt";
Video = false;
Log_file = 'depth_test_2.xlsx';

% Load AUV model parameters
Landmarks = [];
load(pwd + "/Libraries/AUVs/" + model + "/physical.mat");
Sz = 330;

% Read log file data
Test = readmatrix(Log_file);
Test(1,:) = []; %--> Remove header from data
Test(:, 1) = []; %--> Remove numbering column
Test(:, 15) = Test(:, 15) - Test(1, 15); %--> Offset time stamps

% Model setup
dT = 0.1; %--> Initial sampling time for model
n = [Test(1, 1), Test(1, 2), Test(1, 3), Test(1, 4), Test(1, 5), Test(1, 6)]; %--> Initial global pose
v = zeros(1, 6); %--> Initial local velocity
v_dot = zeros(1, 6); %--> Initial acceleration
dT_pwm = round(Test(2, 15) - Test(1, 15), 1); %--> PWM transition sampling time
i = 0; %--> Iteration number

% Create objects
Vehicle = AUV(model, Landmarks, n, v, v_dot); %--> AUV object

% Histories (for records, plots and logs)
NN = [0, n]; %--> Global pose history
VV = [0, v]; %--> Log initial velocity into history vector

% Run model
while(gett(Vehicle) < Test(end, 15))
    % Generate AUV pwm signals "from control node"
    if(gett(Vehicle) >= i * dT_pwm)
        pwm = Test(i + 1, 7 : 14);
        i = i + 1;
    end
                  
    % Update AUV states
    [Vehicle, dT] = Solve(Vehicle, pwm, dT);
                           
    % Update history of states 
    NN = [NN; gett(Vehicle), getn(Vehicle)]; 
    VV = [VV; gett(Vehicle), getv(Vehicle)];
end

% Plots
% Pose (Ground truth and modeled) vs. time
Pose_plot = ["X (m)", "Y (m)", "Z (m)", "Roll (rad)", "Pitch (rad)", "Yaw (rad)"];
for i = 1 : 6
    figure
    plot(NN(:, 1), NN(:, i + 1), 'r', 'LineWidth', 2)
    hold on 
    plot(Test(:, 15), Test(:, i), '--b', 'LineWidth', 2)
    annotation('textbox', [.0 .0 .0 1.0], ...
    'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
    xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
    ylabel(Pose_plot(1, i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
    set (gca, 'fontweight', 'bold', 'FontSize', 18) 
    set(gcf,'units','points','position',[.0, .0, Sz, Sz])
    legend('Model', 'Ground truth', 'Location', 'Best')
end

% Velocity (Ground truth and modeled) vs. time
Vel_plot = ["Surge (m/s)", "Starboard (m/s)", "Heave (m/s)", ...
"About surge (rad/s)", "About starboard (rad/s)", "About heave (rad/s)"];
V_test = [(Test(2 : end, 1) - Test(1 : end - 1, 1)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15)), ...
    (Test(2 : end, 2) - Test(1 : end - 1, 2)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15)), ...
    (Test(2 : end, 3) - Test(1 : end - 1, 3)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15)), ...
    (Test(2 : end, 4) - Test(1 : end - 1, 4)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15)), ...
    (Test(2 : end, 5) - Test(1 : end - 1, 5)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15)), ...
    (Test(2 : end, 6) - Test(1 : end - 1, 6)) ./ (Test(2 : end, 15) - Test(1 : end - 1, 15))]; 
   %--> Rate of change of pose vector
for i = 1 : size(V_test, 1)
    V_test(i, 1 : 6) = (inv(Jacobian(Test(i, 1 : 6))) * V_test(i, 1 : 6).').';
end
for i = 1 : 6
    figure
    plot(VV(:, 1), VV(:, i + 1), 'r', 'LineWidth', 2)
    hold on 
    plot(Test(1 : end - 1, 15), V_test(:, i), '--g', 'LineWidth', 2)
    annotation('textbox', [.0 .0 .0 1.0], ...
    'String', char('a' + i - 1), 'EdgeColor', 'none','FontSize', 18, 'FontWeight', 'bold')
    xlabel('Time (sec.)','FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
    ylabel(Vel_plot(1, i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
    set (gca, 'fontweight', 'bold', 'FontSize', 18) 
    set(gcf,'units','points','position',[.0, .0, Sz, Sz])
    legend('Model', 'Ground truth', 'Location', 'Best')
end

% Trajectory plot
% To be decided

% Additional options
% Run simulink
if(Video == true)
    % Video duration
    dur_vid = gett(Vehicle);

    % State feed to visualization 
    NNx = [NN(:, 1), NN(:, 2)]; NNy = [NN(:, 1), NN(:, 3)];
    NNz = [NN(:, 1), NN(:, 4)]; NNphi = [NN(:, 1), NN(:, 5)];
    NNtht = [NN(:, 1), NN(:, 6)]; NNeps = [NN(:, 1), NN(:, 7)];

    % Run simulation
    sim(pwd + "/Pools/Girls_college_1/test_field.slx");
end

function J = Jacobian(n)
    % Return jacobian relating global state rate of change and local velocity
    J = [cos(n(6))*cos(n(5)) -sin(n(6))*cos(n(4))+cos(n(6))*sin(n(5))*sin(n(4))...
        sin(n(6))*sin(n(4))+cos(n(6))*cos(n(4))*sin(n(5)) 0 0 0;
        sin(n(6))*cos(n(5)) cos(n(6))*cos(n(5))+sin(n(4))*sin(n(5))*sin(n(6)) ...
        -cos(n(6))*sin(n(4))+sin(n(5))*sin(n(6))*cos(n(4)) 0 0 0;
        -sin(n(5)) cos(n(5))*sin(n(4)) cos(n(5))*cos(n(4)) 0 0 0;
        0 0 0 1 sin(n(4))*tan(n(5)) cos(n(4))*tan(n(5))
        0 0 0 0 cos(n(4)) -sin(n(4));
        0 0 0 0 sin(n(4))/cos(n(5)) cos(n(4))/cos(n(5))];
end
