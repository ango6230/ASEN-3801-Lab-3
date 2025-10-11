% Contributors: Owen Storer, Bridger Cushman, Andrew Gonzalez, Sean Sharp
% Assignment: 3801 Lab 3

% House Cleaning
clear;clc;close all

% Toggle Switch
BaseSwitch =	0;
Control1Switch= 1;
Control2Switch= 0;
Control3Switch= 0;
Control4Switch= 0;
Control5Switch= 0;
Control6Switch= 0;
ControlChosenSwitch = 0;
Gyro1Switch =	1; % Gyro1 Auto Frequency = 0.2 [Hz], Current = 0.5 [A]
Gyro2Switch =	0; % Gyro2 Auto Frequency = 0.2 [Hz], Current = 1 [A] Andrew's Preferred Plot
Gyro3Switch =	0; % Gyro3 Auto Frequency = 1 [Hz], Current = 0.5 [A]
Gyro4Switch =	0; % Gyro4 Manual Frequency = 0.2 [Hz], Current = 0.5 [A] Andrew's Preferred Plot
Gyro5Switch =	0; % Gyro5 Manual Frequency = 0.2 [Hz], Current = 1 [A]
Gyro6Switch =	0; % Gyro6 Manual Frequency = 1 [Hz], Current = 0.5 [A]
Rwheel1Switch = 0;
Rwheel2Switch = 0;

%% BASE RUN

% Loading Base Run
[Time_BASE,Data_BASE] = LoadData_BASE('Lab 3 Data/BASE_T10.csv');
p = polyfit(Data_BASE(:,1),Data_BASE(:,2),1);
val_b_BASE = p(2);% This is the bias in the Gyro
val_k_BASE = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_BASE(:,1));

% Finding min/max values for plots and slope
min_BASE = min(Data_BASE(:,1));
max_BASE = max(Data_BASE(:,1));
minPolyline_BASE = min(Polyline);
maxPolyline_BASE = max(Polyline);

Calibrated_Data_BASE = (Data_BASE(:,2) - val_b_BASE)/val_k_BASE;

% Calculating MOI of the Base
torque_const = 25.5/1000; %mNm/Amp
MOI_BASE = (torque_const*mean(Data_BASE(:,3)))/p(1); %mkgm^2


%% CONTROL RUNS

% Loading Control1 Run with k1(proportional) = 50, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL1,Data_CONTROL1] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_50_K2_15_K3_0.csv');
p = polyfit(Data_CONTROL1(:,1),Data_CONTROL1(:,2),1);
val_b_CONTROL1 = p(2);% This is the bias in the Control
val_k_CONTROL1 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_CONTROL1(:,1));

% Finding min/max values for plots and slope
min_CONTROL1 = min(Data_CONTROL1(:,1));
max_CONTROL1 = max(Data_CONTROL1(:,1));
minPolyline_CONTROL1 = min(Polyline);
maxPolyline_CONTROL1 = max(Polyline);

Calibrated_Data_CONTROL1 = (Data_CONTROL1(:,2) - val_b_CONTROL1)/val_k_CONTROL1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control2 Run with k1(proportional) = 100, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL2,Data_CONTROL2] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_15_K3_0.csv');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control3 Run with k1(proportional) = 100, k2(derivative) = 40, k3(integral) = 0
[Time_CONTROL3,Data_CONTROL3] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_40_K3_0.csv');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control4 Run with k1(proportional) = 100, k2(derivative) = -10, k3(integral) = 0
[Time_CONTROL4,Data_CONTROL4] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_neg10_K3_0.csv');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control5 Run with k1(proportional) = 150, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL5,Data_CONTROL5] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_150_K2_15_K3_0.csv');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control6 Run with k1(proportional) = 200, k2(derivative) = 0, k3(integral) = 0
[Time_CONTROL6,Data_CONTROL6] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_200_K2_0_K3_0.csv');

% Loading Control Chosen values of k1 = 3.48 k2 = 16.7
[Time_CONTROL_Chosen, Data_CONTROL_Chosen] = LoadData_CONTROL('Lab 3 Data/10-7-25_CONTROL_kp_3pt48_kd_16pt7.csv');


%% GYRO RUNS

% Loading Gyro1 Run with Auto Frequency = 0.2 [Hz], Current = 0.5 [A]
[Time_GYRO1,Data_GYRO1] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F0pt2_C0pt5.csv');
p = polyfit(Data_GYRO1(:,1),Data_GYRO1(:,2),1);
val_b_GYRO1 = p(2);% This is the bias in the Gyro
val_k_GYRO1 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO1(:,1));

% Finding min/max values for plots and slope
min_GYRO1 = min(Data_GYRO1(:,1));
max_GYRO1 = max(Data_GYRO1(:,1));
minPolyline_GYRO1 = min(Polyline);
maxPolyline_GYRO1 = max(Polyline);

Calibrated_Data_GYRO1 = (Data_GYRO1(:,2) - val_b_GYRO1)/val_k_GYRO1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro2 Run with Auto Frequency = 0.2 [Hz], Current = 1 [A]
[Time_GYRO2,Data_GYRO2] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F0pt2_C1.csv');
p = polyfit(Data_GYRO2(:,1),Data_GYRO2(:,2),1);
val_b_GYRO2 = p(2);% This is the bias in the Gyro
val_k_GYRO2 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO2(:,1));

% Finding min/max values for plots and slope
min_GYRO2 = min(Data_GYRO2(:,1));
max_GYRO2 = max(Data_GYRO2(:,1));
minPolyline_GYRO2 = min(Polyline);
maxPolyline_GYRO2 = max(Polyline);

Calibrated_Data_GYRO2 = (Data_GYRO2(:,2) - val_b_GYRO2)/val_k_GYRO2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro3 Run with Auto Frequency = 1 [Hz], Current = 0.5 [A]
[Time_GYRO3,Data_GYRO3] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F1_C0pt5.csv');
p = polyfit(Data_GYRO3(:,1),Data_GYRO3(:,2),1);
val_b_GYRO3 = p(2);% This is the bias in the Gyro
val_k_GYRO3 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO3(:,1));

% Finding min/max values for plots and slope
min_GYRO3 = min(Data_GYRO3(:,1));
max_GYRO3 = max(Data_GYRO3(:,1));
minPolyline_GYRO3 = min(Polyline);
maxPolyline_GYRO3 = max(Polyline);

Calibrated_Data_GYRO3 = (Data_GYRO3(:,2) - val_b_GYRO3)/val_k_GYRO3;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro4 Run with Manual Frequency = 0.2 [Hz], Current = 0.5 [A]
[Time_GYRO4,Data_GYRO4] = LoadData_GYRO('Lab 3 Data/GYRO_MAN_F0pt2_C0pt5.csv');
p = polyfit(Data_GYRO4(:,1),Data_GYRO4(:,2),1);
val_b_GYRO4 = p(2);% This is the bias in the Gyro
val_k_GYRO4 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO4(:,1));

% Finding min/max values for plots and slope
min_GYRO4 = min(Data_GYRO4(:,1));
max_GYRO4 = max(Data_GYRO4(:,1));
minPolyline_GYRO4 = min(Polyline);
maxPolyline_GYRO4 = max(Polyline);

Calibrated_Data_GYRO4 = (Data_GYRO4(:,2) - val_b_GYRO4)/val_k_GYRO4;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro5 Run with Manual Frequency = 0.2 [Hz], Current = 1 [A]
[Time_GYRO5,Data_GYRO5] = LoadData_GYRO('Lab 3 Data/GYRO_MAN_F0pt2_C1.csv');
p = polyfit(Data_GYRO5(:,1),Data_GYRO5(:,2),1);
val_b_GYRO5 = p(2);% This is the bias in the Gyro
val_k_GYRO5 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO5(:,1));

% Finding min/max values for plots and slope
min_GYRO5 = min(Data_GYRO5(:,1));
max_GYRO5 = max(Data_GYRO5(:,1));
minPolyline_GYRO5 = min(Polyline);
maxPolyline_GYRO5 = max(Polyline);

Calibrated_Data_GYRO5 = (Data_GYRO5(:,2) - val_b_GYRO5)/val_k_GYRO5;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro6 Run with Manual Frequency = 1 [Hz], Current = 0.5 [A]
[Time_GYRO6,Data_GYRO6] = LoadData_GYRO('Lab 3 Data/GYRO_MAN_F1_C0pt5.csv');
p = polyfit(Data_GYRO6(:,1),Data_GYRO6(:,2),1);
val_b_GYRO6 = p(2);% This is the bias in the Gyro
val_k_GYRO6 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO5(:,1));

% Finding min/max values for plots and slope
min_GYRO6 = min(Data_GYRO6(:,1));
max_GYRO6 = max(Data_GYRO6(:,1));
minPolyline_GYRO6 = min(Polyline);
maxPolyline_GYRO6 = max(Polyline);

%%
angular_rate_error_1 = (abs(Data_GYRO1(:,2)) - abs(Calibrated_Data_GYRO1));
angular_rate_error_3 = (abs(Data_GYRO3(:,2)) - abs(Calibrated_Data_GYRO3));

figure();
plot(Time_GYRO1,angular_rate_error_1);
hold on;
xlabel('Time(s)');
ylabel('Angular Rate Measurement Error');
title('Time History Of Angular Rate Measurement Error for Data Set 1');

figure();
plot(Time_GYRO3,angular_rate_error_3);
hold on;
xlabel('Time(s)');
ylabel('Angular Rate Measurement Error');
title('Time History Of Angular Rate Measurement Error for Data Set 3');

%true angular position
theta0_1 = 0;
theta_1= theta0_1 + cumtrapz(Time_GYRO1, Data_GYRO1(:,2));
theta0_3 = 0;
theta_3= theta0_3 + cumtrapz(Time_GYRO3, Data_GYRO3(:,2));

%measured angular position
phi0_1 = 0;
phi_1 = phi0_1 + cumtrapz(Time_GYRO1, Calibrated_Data_GYRO1(:,1));
phi0_3 = 0;
phi_3 = phi0_3 + cumtrapz(Time_GYRO3, Calibrated_Data_GYRO3(:,1));


figure();
plot(Time_GYRO1, theta_1);
hold on;
plot(Time_GYRO1, -phi_1);
ylabel('Rate Measurement'); %what units 
xlabel('Time(s)');
title('Time History of True Angular Position(Encoder) and Measured Angular Position(Gyro) For Data Set 1')
legend('True Angular Position','Measured Angular Position');

figure();
plot(Time_GYRO3, theta_3);
hold on;
plot(Time_GYRO3, -phi_3);
ylabel('Rate Measurement'); %what units 
xlabel('Time(s)');
title('Time History of True Angular Position(Encoder) and Measured Angular Position(Gyro) For Data Set 3')
legend('True Angular Position','Measured Angular Position');



Calibrated_Data_GYRO6 = (Data_GYRO6(:,2) - val_b_GYRO6)/val_k_GYRO6;

%% REACTION WHEEL RUNS

% Loading Reaction Wheel Test Run with Torque = 0.2
[Time_RWHEEL1,Data_RWHEEL1] = LoadData_RWHEEL('Lab 3 Data/RWHEEL_T0pt2_Test.csv');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Reaction Wheel Run with Torque = 20
[Time_RWHEEL2,Data_RWHEEL2] = LoadData_RWHEEL('Lab 3 Data/RWHEEL_T20.csv');

%% Plotting Section

%{
figure()
hold on
scatter(Data_BASE(:,1),Data_BASE(:,2),'.');
grid on
plot([min_BASE,max_BASE],[maxPolyline_BASE,minPolyline_BASE],'--')
yline(val_b_BASE,'--');
hold off

figure()
hold on
plot(Time_BASE,Data_BASE(:,1))
plot(Time_BASE,Calibrated_Data_BASE)
hold off
%}

if Control1Switch == 1
    figure(name='Control Kp=50 Kd=15')
    axis square
    hold on
    grid on
    plot(Time_CONTROL1-Time_CONTROL1(1),Calibrated_Data_CONTROL1)
    plot(Time_CONTROL1-Time_CONTROL1(1), Data_CONTROL1(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Control2Switch == 1
    figure(name='Control Kp=100 Kd=15')
    axis square
    hold on
    grid on
    plot(Time_CONTROL2-Time_CONTROL2(1), unwrap(Data_CONTROL2(:,2))-Data_CONTROL2(1,2))
    plot(Time_CONTROL2-Time_CONTROL2(1), Data_CONTROL2(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Control3Switch == 1
    figure(name='Control Kp=100 Kd=40')
    axis square
    hold on
    grid on
    plot(Time_CONTROL3-Time_CONTROL3(1), unwrap(Data_CONTROL3(:,2))-Data_CONTROL3(1,2))
    plot(Time_CONTROL3-Time_CONTROL3(1), Data_CONTROL3(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Control4Switch == 1
    figure(name='Control Kp=100 Kd=-10')
    axis square
    hold on
    grid on
    plot(Time_CONTROL4-Time_CONTROL4(1), unwrap(Data_CONTROL4(:,2))-Data_CONTROL4(1,2))
    plot(Time_CONTROL4-Time_CONTROL4(1), Data_CONTROL4(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Control5Switch == 1
    figure(name='Control Kp=150 Kd=15')
    axis square
    hold on
    grid on
    plot(Time_CONTROL5-Time_CONTROL5(1), unwrap(Data_CONTROL5(:,2))-Data_CONTROL5(1,2))
    plot(Time_CONTROL5-Time_CONTROL5(1), Data_CONTROL5(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Control6Switch == 1
    figure(name='Control Kp=200 Kd=0')
    hold on
    grid on
    plot(Time_CONTROL6-Time_CONTROL6(1), unwrap(Data_CONTROL6(:,2)-Data_CONTROL6(1,2)))
    plot(Time_CONTROL6-Time_CONTROL6(1), Data_CONTROL6(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if ControlChosenSwitch == 1
    figure(name='Control Kp=3.48 Kd=16.7')
    hold on
    grid on
    plot(Time_CONTROL_Chosen-Time_CONTROL_Chosen(1), unwrap(Data_CONTROL_Chosen(:,2))-Data_CONTROL_Chosen(1,2))
    plot(Time_CONTROL_Chosen-Time_CONTROL_Chosen(1), Data_CONTROL_Chosen(:,1))
    xlabel('Time [s]')
    ylabel('Position [rad]')
    hold off
end
if Gyro1Switch == 1
	figure(name='Gyro 1 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO1(:,2),Data_GYRO1(:,1),'.');
	plot([maxPolyline_GYRO1,minPolyline_GYRO1],[min_GYRO1,max_GYRO1],'--')
	yline(val_b_GYRO1,'--');
	hold off

	figure(name='Gyro 1 Time Plot')
	hold on
	plot(Time_GYRO1,Data_GYRO1(:,1))
	plot(Time_GYRO1,Calibrated_Data_GYRO1)
	%plot(Time_GYRO1,Data_GYRO1)
	hold off
end

if Gyro2Switch == 1
	% Figure plot for Problem 3.1 part c-1
	figure(name='Raw Unflipped Gyro 2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,1),'--',LineWidth=2)
	plot(Time_GYRO2,Data_GYRO2(:,2),'r',LineWidth=0.8)
	title('Raw Gyro Data Overlaid with Encoder Data Before Calibration',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Raw Gyro Data')
	hold off
	print('Problem_3_1c1-1','-dpng','-r600')

	% Figure plot for Problem 3.1 part c-1
	figure(name='Raw Flipped Gyro 2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,1),'--',LineWidth=2)
	plot(Time_GYRO2,-1 .* Data_GYRO2(:,2),'r',LineWidth=0.8)
	title('Raw Gyro Data Overlaid with Encoder Data Before Calibration',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Raw Gyro Data')
	hold off
	print('Problem_3_1c1-2','-dpng','-r600')

	% Figure plot for Problem 3.1 part c-2
	figure(name='Gyro 2 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO2(:,2),Data_GYRO2(:,1),'.');
	plot([maxPolyline_GYRO2,minPolyline_GYRO2],[min_GYRO2,max_GYRO2],'--r',LineWidth=2)
	yline(val_b_GYRO2,'--');
	title('Gyro Calibration Data',FontSize=18)
	xlabel('Encoder Rate Measurement [rad/s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Data','Adjusted Scale Factor \itK','Bias \itb')
	hold off
	print('Problem_3_1c2-1','-dpng','-r600')
	
	% Figure plot for Problem 3.1 part c-2
	figure(name='Calibrated Gyro 2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,1),'--',LineWidth=2)
	plot(Time_GYRO2,Calibrated_Data_GYRO2,'r',LineWidth=0.8)
	title('Calibrated Gyro Data Overlaid with Encoder Data',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data')
	hold off
	print('Problem_3_1c2-2','-dpng','-r600')
end
	
if Gyro3Switch == 1
	figure(name='Gyro 3 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO3(:,2),Data_GYRO3(:,1),'.');
	plot([maxPolyline_GYRO3,minPolyline_GYRO3],[min_GYRO3,max_GYRO3],'--')
	yline(val_b_GYRO3,'--');
	hold off
	
	figure(name='Gyro 3 Time Plot')
	hold on
	plot(Time_GYRO3,Data_GYRO3(:,1))
	plot(Time_GYRO3,Calibrated_Data_GYRO3)
	%plot(Time_GYRO3,Data_GYRO3)
	hold off
end
	
if Gyro4Switch == 1
	figure(name='Gyro 4 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO4(:,2),Data_GYRO4(:,1),'.');
	plot([maxPolyline_GYRO4,minPolyline_GYRO4],[min_GYRO4,max_GYRO4],'--')
	yline(val_b_GYRO4,'--');
	hold off
	
	figure(name='Gyro 4 Time Plot')
	hold on
	plot(Time_GYRO4,Data_GYRO4(:,1))
	plot(Time_GYRO4,Calibrated_Data_GYRO4)
	%plot(Time_GYRO4,Data_GYRO4)
	hold off
end
	
if Gyro5Switch == 1
	figure(name='Gyro 5 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO5(:,2),Data_GYRO5(:,1),'.');
	plot([maxPolyline_GYRO5,minPolyline_GYRO5],[min_GYRO5,max_GYRO5],'--')
	yline(val_b_GYRO5,'--');
	hold off
	
	figure(name='Gyro 5 Time Plot')
	hold on
	plot(Time_GYRO5,Data_GYRO5(:,1))
	plot(Time_GYRO5,Calibrated_Data_GYRO5)
	%plot(Time_GYRO5,Data_GYRO5)
	hold off
end
	
if Gyro6Switch == 1
	figure(name='Gyro 6 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO6(:,2),Data_GYRO6(:,1),'.');
	plot([maxPolyline_GYRO6,minPolyline_GYRO6],[min_GYRO6,max_GYRO6],'--')
	yline(val_b_GYRO6,'--');
	hold off
	
	figure(name='Gyro 6 Time Plot')
	hold on
	plot(Time_GYRO6,Data_GYRO6(:,1))
	plot(Time_GYRO6,Calibrated_Data_GYRO6)
	%plot(Time_GYRO6,Data_GYRO6)
	hold off
end

%% Functions Section

function [Time,Data] = LoadData_BASE(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data(:,1) = data(:,1)/1000; % Convert from milliseconds to seconds
data(:,3) = data(:,3)/60 * 2 * pi; % Convert from RPM to rad/s
data = unique(data, 'rows'); % Remove duplicate entries
Time = data(:, 1); % Extracting the time column
Data = data(:, 2:end); % Extracting the data columns

end

function [Time,Data] = LoadData_CONTROL(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data(:,1) = data(:,1)/1000; % Convert from milliseconds to seconds
data(1,:) = [];
data = unique(data, 'rows'); % Remove duplicate entries
Time = data(:, 1); % Extracting the time column
Data = data(:, 2:end); % Extracting the data columns

end

function [Time,Data] = LoadData_GYRO(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data(:,3) = data(:,3)/60 * 2 * pi; % Convert from RPM to rad/s
truncate = find(data(:,2) >= 8.7,1);
data(1:truncate,:) = [];
mask = data(:,2) >= 8.7;
data(mask,:) = [];
data(1,:) = [];
data = unique(data, 'rows'); % Remove duplicate entries
Time = data(:, 1); % Extracting the time column
Time(:,1) = Time(:,1) - Time(1);
Data = data(:, 2:end); % Extracting the data columns

end

function [Time,Data] = LoadData_RWHEEL(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data(:,1) = data(:,1)/1000; % Convert from milliseconds to seconds
data(1,:) = [];
data = unique(data, 'rows'); % Remove duplicate entries
Time = data(:, 1); % Extracting the time column
Time(:,1) = Time(:,1) - Time(1);
Data = data(:, 2:end); % Extracting the data columns

end
