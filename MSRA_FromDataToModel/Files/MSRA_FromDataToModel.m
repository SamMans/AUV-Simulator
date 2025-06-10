% This script is the programmatic workflow associated with generating the
% model associated with the Blue Robotics T200 Thruster as shown in the
% Robotics Arena Video MATLAB and Simulink Robotics Arena : From Data To
% Model.

% Copyright 2018 The MathWorks, Inc.

%% Import Data into the workspace and save it to the iddata object 
filenameSine = 'T200_Sine_0-10_Hz_1600-1900_us.csv';
filenameSquare = 'T200_Square_0-10_Hz_1600-1900_us.csv';
startRow = 2;
endRow = 15002;
Ts = 0.002; %Data Collection sample time
[input_sn,force_sn] = importData(filenameSine, startRow, endRow);
[input_sq,force_sq] = importData(filenameSquare, startRow, endRow);
sine = iddata(force_sn,input_sn,Ts);
square = iddata(force_sq,input_sq,Ts);

clearvars filenameSine filenameSquare startRow endRow

figure 
subplot(2,1,1)
plot(sine)
subplot(2,1,2)
plot(square)

%% Pre-processing

% Remove Means 
sined = detrend(sine,0);
squared = detrend(square,0);

% Filter 
sinedf = idfilt(sined,5,0.07394);
squaredf = idfilt(squared,5,0.064399);
%% Identify ss1 Model from the video 

Options = n4sidOptions;
Options.Display = 'on';
Options.N4Weight = 'CVA';
Options.N4Horizon = 'Auto';

ss1 = n4sid(squared,1:10,'Ts',0.002,Options);

%Validation
figure
subplot(3,1,1)
compare(sined,ss1);

%% Identify ss2 Model from the video

Options = n4sidOptions;
Options.Display = 'on';
Options.EnforceStability = 1;
Options.N4Weight = 'CVA';
Options.N4Horizon = 'Auto';

ss2 = n4sid(squaredf,1:10,'Ts',0.002,Options);

%Validation
subplot(3,1,2)
compare(sinedf,ss2);

%% Identify Transfer function model from the video 

Options = tfestOptions;
Options.Display = 'on';
Options.WeightingFilter = [];

tf1 = tfest(squaredf,4,2,Options);

% Validation
subplot(3,1,3)
compare(sinedf,tf1);




