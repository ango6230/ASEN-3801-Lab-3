% Contributors: Owen Storer, Bridger Cushman, Andrew Gonzalez, Sean Sharp
% Assignment: 3801 Lab 3

% Housekeeping
clear;clc;close all

%% Loading all unique data files in

% Loading Base Run
[Time_BASE,Data_BASE] = LoadData_BASE('Lab 3 Data/BASE_T10.csv');
%time = time/1000; % Time array now in seconds

% Loading Control1 Run with k1 = 50, k2 = 15, k3 = 0
[Time_CONTROL1,Data_CONTROL1] = LoadData('Lab 3 Data/CONTROL_K1_50_K2_15_K3_0.csv');

% Loading Control2 Run with k1 = 100, k2 = 15, k3 = 0
[Time_CONTROL2,Data_CONTROL2] = LoadData('Lab 3 Data/CONTROL_K1_100_K2_15_K3_0.csv');

% Loading Control3 Run with k1 = 100, k2 = 40, k3 = 0
[Time_CONTROL3,Data_CONTROL3] = LoadData('Lab 3 Data/CONTROL_K1_100_K2_40_K3_0.csv');

% Loading Control4 Run with k1 = 100, k2 = -10, k3 = 0
[Time_CONTROL4,Data_CONTROL4] = LoadData('Lab 3 Data/CONTROL_K1_100_K2_neg10_K3_0.csv');

% Loading Control5 Run with k1 = 150, k2 = 15, k3 = 0
[Time_CONTROL5,Data_CONTROL5] = LoadData('Lab 3 Data/CONTROL_K1_150_K2_15_K3_0.csv');

% Loading Control6 Run with k1 = 200, k2 = 0, k3 = 0
[Time_CONTROL6,Data_CONTROL6] = LoadData('Lab 3 Data/CONTROL_K1_200_K2_0_K3_0.csv');

% Loading Gyro1 Run with Auto F = 0.2, C = 0.5
[Time_GYRO1,Data_GYRO1] = LoadData('Lab 3 Data/GYRO_AUTO_F0pt2_C0pt5.csv');

% Loading Gyro2 Run with Auto F = 0.2, C = 1
[Time_GYRO2,Data_GYRO2] = LoadData('Lab 3 Data/GYRO_AUTO_F0pt2_C1.csv');

% Loading Gyro3 Run with Auto F = 1, C = 0.5
[Time_GYRO3,Data_GYRO3] = LoadData('Lab 3 Data/GYRO_AUTO_F1_C0pt5.csv');

% Loading Gyro4 Run with Manual F = 0.2, C = 0.5
[Time_GYRO4,Data_GYRO4] = LoadData('Lab 3 Data/GYRO_MAN_F0pt2_C0pt5.csv');

% Loading Gyro5 Run with Manual F = 0.2, C = 1
[Time_GYRO5,Data_GYRO5] = LoadData('Lab 3 Data/GYRO_MAN_F0pt2_C1.csv');

% Loading Gyro6 Run with Manual F = 1, C = 0.5
[Time_GYRO6,Data_GYRO6] = LoadData('Lab 3 Data/GYRO_MAN_F1_C0pt5.csv');

% Loading Reaction Wheel Test Run with T = 0.2
[Time_RWHEEL1,Data_RWHEEL1] = LoadData('Lab 3 Data/RWHEEL_T0pt2_Test.csv');

% Loading Reaction Wheel Run with T = 20
[Time_RWHEEL2,Data_RWHEEL2] = LoadData('Lab 3 Data/RWHEEL_T20.csv');

%% Functions Section

function [Time,Data] = LoadData_BASE(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data = unique(data, 'rows');
Time = data(:, 1); % Extracting the time column
Time = Time/1000; % Converting from milliseconds to seconds
Data = data(:, 2:end); % Extracting the data columns

end

function [Time,Data] = LoadData(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data = unique(data, 'rows');
Time = data(:, 1); % Extracting the time column
Time = Time(2:end ,1) - Time(2);
Time = Time/1000; % Converting from milliseconds to seconds
Data = data(:, 2:end); % Extracting the data columns

end