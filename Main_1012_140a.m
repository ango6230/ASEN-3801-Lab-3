% Contributors: Owen Storer, Bridger Cushman, Andrew Gonzalez, Sean Sharp
% Assignment: 3801 Lab 3

% House Cleaning
clear;clc;close all

% Toggle Switch 0=off, 1=on
BaseSwitch =	0;
Control1Switch= 0;
Control2Switch= 0;
Control3Switch= 0;
Control4Switch= 0;
Control5Switch= 0;
Control6Switch= 0;
ControlChosenSwitch = 0;
Gyro1Switch =	0; % Gyro1 Auto Frequency = 0.2 [Hz], Current = 0.5 [A]
Gyro2Switch =	0; % Gyro2 Auto Frequency = 0.2 [Hz], Current = 1 [A] Andrew's Preferred Plot
Gyro3Switch =	0; % Gyro3 Auto Frequency = 1 [Hz], Current = 0.5 [A]
Gyro4Switch =	0; % Gyro4 Manual Frequency = 0.2 [Hz], Current = 0.5 [A]
Gyro5Switch =	0; % Gyro5 Manual Frequency = 0.2 [Hz], Current = 1 [A]
Gyro6Switch =	0; % Gyro6 Manual Frequency = 1 [Hz], Current = 0.5 [A] Andrew's Preferred Plot
Rwheel1Switch = 0; %T=5mNm
Rwheel2Switch = 0; %T=10mNm
Rwheel3Switch = 0; %T=15mNm
Rwheel4Switch = 0; %T=20mNm

%% BASE RUN

TorqueConst_BASE = 25.5; % [mNm/A]

% Loading Base Run
[Time_BASE,Data_BASE] = LoadData_BASE('Lab 3 Data/BASE_T10.csv');
% Converting Amps to Torque [Nm]
Torque_BASE = (Data_BASE(:,3) * TorqueConst_BASE)/1000;
%Torque_BASE = nonzeros(Torque_BASE);
Torque_BASE = mean(Torque_BASE);
% Converting angular velocity to angular acceleration [rad/s^2]
alpha_BASE = gradient((pi/180)*Data_BASE(:,2), Time_BASE);
%alpha_BASE = nonzeros(alpha_BASE);
alpha_BASE = mean(alpha_BASE);
% Calculating MOI of the Base [kg*m^2]
MOI_BASE = Torque_BASE/alpha_BASE;

%{
%p = polyfit(Data_BASE(:,1),Data_BASE(:,2),1);
p = polyfit(Time_BASE,Data_BASE(:,2),1);
val_b_BASE = p(2);% This is the bias in the Gyro
val_k_BASE = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Time_BASE);

plot(Time_BASE,Polyline)

% Finding min/max values for plots and slope
min_BASE = min(Data_BASE(:,1));
max_BASE = max(Data_BASE(:,1));
minPolyline_BASE = min(Polyline);
maxPolyline_BASE = max(Polyline);

% Calculating MOI of the Base [kg*m^2]
%MOI_BASE = (TorqueConst_BASE * mean(Data_BASE(:,3)))/mean(alpha_BASE);
MOI_BASE = Torque_BASE/alpha_BASE;
%}

%% CONTROL RUNS

% Loading Control1 Run with k1(proportional) = 50, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL1,Data_CONTROL1] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_50_K2_15_K3_0.csv');
Time_CONTROL1 = Time_CONTROL1 - Time_CONTROL1(1);
unwrapped_CONTROL1 = unwrap(Data_CONTROL1(:,2));
%shift_CONTROL1 = unwrapped_CONTROL1 - unwrapped_CONTROL1(end) + Data_CONTROL1(end,1);
shift_CONTROL1 = unwrapped_CONTROL1 - unwrapped_CONTROL1(end) + 0.25;
wrapped_CONTROL1 = wrapToPi(shift_CONTROL1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control2 Run with k1(proportional) = 100, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL2,Data_CONTROL2] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_15_K3_0.csv');
Time_CONTROL2 = Time_CONTROL2 - Time_CONTROL2(1);
unwrapped_CONTROL2 = unwrap(Data_CONTROL2(:,2));
%shift_CONTROL2 = unwrapped_CONTROL2 - unwrapped_CONTROL2(end) + Data_CONTROL2(end,1);
shift_CONTROL2 = unwrapped_CONTROL2 - unwrapped_CONTROL2(end) + 0.25;
wrapped_CONTROL2 = wrapToPi(shift_CONTROL2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control3 Run with k1(proportional) = 100, k2(derivative) = 40, k3(integral) = 0
[Time_CONTROL3,Data_CONTROL3] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_40_K3_0.csv');
Time_CONTROL3 = Time_CONTROL3 - Time_CONTROL3(1);
unwrapped_CONTROL3 = unwrap(Data_CONTROL3(:,2));
%shift_CONTROL3 = unwrapped_CONTROL3 - unwrapped_CONTROL3(end) + Data_CONTROL3(end,1);
shift_CONTROL3 = unwrapped_CONTROL3 - unwrapped_CONTROL3(end) + 0.25;
wrapped_CONTROL3 = wrapToPi(shift_CONTROL3);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control4 Run with k1(proportional) = 100, k2(derivative) = -10, k3(integral) = 0
[Time_CONTROL4,Data_CONTROL4] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_100_K2_neg10_K3_0.csv');
Time_CONTROL4 = Time_CONTROL4 - Time_CONTROL4(1);
unwrapped_CONTROL4 = unwrap(Data_CONTROL4(:,2));
%shift_CONTROL4 = unwrapped_CONTROL4 - unwrapped_CONTROL4(end) + Data_CONTROL4(end,1);
shift_CONTROL4 = unwrapped_CONTROL4 - unwrapped_CONTROL4(end) + 0.25;
wrapped_CONTROL4 = wrapToPi(shift_CONTROL4);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control5 Run with k1(proportional) = 150, k2(derivative) = 15, k3(integral) = 0
[Time_CONTROL5,Data_CONTROL5] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_150_K2_15_K3_0.csv');
Time_CONTROL5 = Time_CONTROL5 - Time_CONTROL5(1);
unwrapped_CONTROL5 = unwrap(Data_CONTROL5(:,2));
%shift_CONTROL5 = unwrapped_CONTROL5 - unwrapped_CONTROL5(end) + Data_CONTROL5(end,1);
shift_CONTROL5 = unwrapped_CONTROL5 - unwrapped_CONTROL5(end) + 0.25;
wrapped_CONTROL5 = wrapToPi(shift_CONTROL5);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control6 Run with k1(proportional) = 200, k2(derivative) = 0, k3(integral) = 0
[Time_CONTROL6,Data_CONTROL6] = LoadData_CONTROL('Lab 3 Data/CONTROL_K1_200_K2_0_K3_0.csv');
Time_CONTROL6 = Time_CONTROL6 - Time_CONTROL6(1);
unwrapped_CONTROL6 = unwrap(Data_CONTROL6(:,2));
%shift_CONTROL6 = unwrapped_CONTROL6 - unwrapped_CONTROL6(end) + Data_CONTROL6(end,1);
shift_CONTROL6 = unwrapped_CONTROL6 - unwrapped_CONTROL6(end);
wrapped_CONTROL6 = wrapToPi(shift_CONTROL6);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Control Chosen values of k1 = 3.48 k2 = 16.7
[Time_CONTROL_Chosen, Data_CONTROL_Chosen] = LoadData_CONTROL('Lab 3 Data/10-7-25_CONTROL_kp_3pt48_kd_16pt6.csv');
Time_CONTROL_Chosen = Time_CONTROL_Chosen - Time_CONTROL_Chosen(1);

%% GYRO RUNS

% Loading Gyro1 Run with Auto Frequency = 0.2 [Hz], Current = 0.5 [A]
[Time_GYRO1,Data_GYRO1] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F0pt2_C0pt5.csv');
p = polyfit(Data_GYRO1(:,2),Data_GYRO1(:,1),1);
val_b_GYRO1 = p(2);% This is the bias in the Gyro
val_k_GYRO1 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO1(:,2));

% Finding min/max values for plots and slope
min_GYRO1 = min(Data_GYRO1(:,2));
max_GYRO1 = max(Data_GYRO1(:,2));
minPolyline_GYRO1 = min(Polyline);
maxPolyline_GYRO1 = max(Polyline);

Calibrated_Data_GYRO1 = (Data_GYRO1(:,1) - val_b_GYRO1)/val_k_GYRO1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro2 Run with Auto Frequency = 0.2 [Hz], Current = 1 [A]
[Time_GYRO2,Data_GYRO2] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F0pt2_C1.csv');
p = polyfit(Data_GYRO2(:,2),Data_GYRO2(:,1),1);
val_b_GYRO2 = p(2);% This is the bias in the Gyro
val_k_GYRO2 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO2(:,2));

% Finding min/max values for plots and slope
min_GYRO2 = min(Data_GYRO2(:,2));
max_GYRO2 = max(Data_GYRO2(:,2));
minPolyline_GYRO2 = min(Polyline);
maxPolyline_GYRO2 = max(Polyline);

Calibrated_Data_GYRO2 = (Data_GYRO2(:,1) - val_b_GYRO2)/val_k_GYRO2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro3 Run with Auto Frequency = 1 [Hz], Current = 0.5 [A]
[Time_GYRO3,Data_GYRO3] = LoadData_GYRO('Lab 3 Data/GYRO_AUTO_F1_C0pt5.csv');
p = polyfit(Data_GYRO3(:,2),Data_GYRO3(:,1),1);
val_b_GYRO3 = p(2);% This is the bias in the Gyro
val_k_GYRO3 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO3(:,2));

% Finding min/max values for plots and slope
min_GYRO3 = min(Data_GYRO3(:,2));
max_GYRO3 = max(Data_GYRO3(:,2));
minPolyline_GYRO3 = min(Polyline);
maxPolyline_GYRO3 = max(Polyline);

Calibrated_Data_GYRO3 = (Data_GYRO3(:,1) - val_b_GYRO3)/val_k_GYRO3;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro4 Run with Manual Frequency = 0.2 [Hz], Current = 0.5 [A]
[Time_GYRO4,Data_GYRO4] = LoadData_GYRO('Lab 3 Data/10-7-25_GYRO_MAN_F0pt2_C0pt5.csv');
p = polyfit(Data_GYRO4(:,2),Data_GYRO4(:,1),1);
val_b_GYRO4 = p(2);% This is the bias in the Gyro
val_k_GYRO4 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO4(:,2));

% Finding min/max values for plots and slope
min_GYRO4 = min(Data_GYRO4(:,2));
max_GYRO4 = max(Data_GYRO4(:,2));
minPolyline_GYRO4 = min(Polyline);
maxPolyline_GYRO4 = max(Polyline);

Calibrated_Data_GYRO4 = (Data_GYRO4(:,1) - val_b_GYRO4)/val_k_GYRO4;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro5 Run with Manual Frequency = 0.2 [Hz], Current = 1 [A]
[Time_GYRO5,Data_GYRO5] = LoadData_GYRO('Lab 3 Data/10-7-25_GYRO_MAN_F0pt2_C1.csv');
p = polyfit(Data_GYRO5(:,2),Data_GYRO5(:,1),1);
val_b_GYRO5 = p(2);% This is the bias in the Gyro
val_k_GYRO5 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO5(:,2));

% Finding min/max values for plots and slope
min_GYRO5 = min(Data_GYRO5(:,2));
max_GYRO5 = max(Data_GYRO5(:,2));
minPolyline_GYRO5 = min(Polyline);
maxPolyline_GYRO5 = max(Polyline);

Calibrated_Data_GYRO5 = (Data_GYRO5(:,1) - val_b_GYRO5)/val_k_GYRO5;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Gyro6 Run with Manual Frequency = 1 [Hz], Current = 0.5 [A]
[Time_GYRO6,Data_GYRO6] = LoadData_GYRO('Lab 3 Data/10-7-25_GYRO_MAN_F1_C0pt5.csv');
p = polyfit(Data_GYRO6(:,2),Data_GYRO6(:,1),1);
val_b_GYRO6 = p(2);% This is the bias in the Gyro
val_k_GYRO6 = p(1);% This is the scale factor for calibration
Polyline = polyval(p,Data_GYRO5(:,2));

% Finding min/max values for plots and slope
min_GYRO6 = min(Data_GYRO6(:,2));
max_GYRO6 = max(Data_GYRO6(:,2));
minPolyline_GYRO6 = min(Polyline);
maxPolyline_GYRO6 = max(Polyline);

Calibrated_Data_GYRO6 = (Data_GYRO6(:,1) - val_b_GYRO6)/val_k_GYRO6;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculating mean and std deviation of b val and k val

val_b_array = [val_b_GYRO1, val_b_GYRO2, val_b_GYRO3, val_b_GYRO4,...
	val_b_GYRO5,val_b_GYRO6];
val_k_array = [val_k_GYRO1, val_k_GYRO2, val_k_GYRO3, val_k_GYRO4,...
	val_k_GYRO5,val_k_GYRO6];

val_b_avg = mean(val_b_array);
val_k_avg = mean(val_k_array);

val_b_std = std(val_b_array);
val_k_std = std(val_k_array);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculating Angular Rate Error (3.1_c4ii)

%angular_rate_error_GYRO1 = abs(Calibrated_Data_GYRO1) - abs(Data_GYRO1(:,2));
%angular_rate_error_GYRO3 = abs(Calibrated_Data_GYRO3) - abs(Data_GYRO3(:,2));
%angular_rate_error_GYRO1 = Calibrated_Data_GYRO1 - Data_GYRO1(:,2);
%angular_rate_error_GYRO3 = Calibrated_Data_GYRO3 - Data_GYRO3(:,2);

% Calculating Angular Position Error (3.1_c4iv)
%angular_position_error_GYRO1 = abs(theta_GYRO1) - abs(phi_GYRO1);
%angular_position_error_GYRO3 = abs(theta_GYRO3) - abs(phi_GYRO3);
%angular_position_error_GYRO1 = sqrt(theta_GYRO1.^2 + phi_GYRO1.^2);
%angular_position_error_GYRO3 = sqrt(theta_GYRO3.^2 + phi_GYRO3.^2);

enc_u  = unwrap(Data_GYRO1(:,2));
gyro_u = unwrap(Calibrated_Data_GYRO1);
enc0   = enc_u  - enc_u(1);
gyro0  = gyro_u - gyro_u(1);
angular_rate_error_GYRO1 = wrapToPi(gyro0 - enc0);

enc_u  = unwrap(Data_GYRO3(:,2));
gyro_u = unwrap(Calibrated_Data_GYRO3);
enc0   = enc_u  - enc_u(1);
gyro0  = gyro_u - gyro_u(1);
angular_rate_error_GYRO3 = wrapToPi(gyro0 - enc0);

% True Angular Position (3.1_c4iii)
theta0_GYRO1 = 0;
theta_GYRO1 = theta0_GYRO1 + cumtrapz(Time_GYRO1, Data_GYRO1(:,2));
theta0_GYRO3 = 0;
theta_GYRO3 = theta0_GYRO3 + cumtrapz(Time_GYRO3, Data_GYRO3(:,2));

% Measured Angular Position (3.1_c4iii)
phi0_GYRO1 = 0;
phi_GYRO1 = phi0_GYRO1 + cumtrapz(Time_GYRO1, Calibrated_Data_GYRO1);
phi0_GYRO3 = 0;
phi_GYRO3 = phi0_GYRO3 + cumtrapz(Time_GYRO3, Calibrated_Data_GYRO3);

theta_enc_u  = unwrap(theta_GYRO1);
theta_gyro_u = unwrap(phi_GYRO1);
theta_enc0   = theta_enc_u  - theta_enc_u(1);
theta_gyro0  = theta_gyro_u - theta_gyro_u(1);
angular_position_error_GYRO1 = wrapToPi(theta_gyro0 - theta_enc0);

theta_enc_u  = unwrap(theta_GYRO3);
theta_gyro_u = unwrap(phi_GYRO3);
theta_enc0   = theta_enc_u  - theta_enc_u(1);
theta_gyro0  = theta_gyro_u - theta_gyro_u(1);
angular_position_error_GYRO3 = wrapToPi(theta_gyro0 - theta_enc0);


%% REACTION WHEEL RUNS

torque_const_RW = 33.5;

% Loading Reaction Wheel Test Run with Torque = 5
[Time_RW1,Data_RW1] = LoadData_RWHEEL('Lab 3 Data/10-7-25_RWHEEL_T5.csv');
p = polyfit(Time_RW1,Data_RW1(:,2),1);
val_b_RW1 = p(2);% This is the bias in the RW
val_k_RW1 = p(1);% This is the scale factor for calibration
Polyline_RW1 = polyval(p,Time_RW1);

Calibrated_Data_RW1 = (Data_RW1(:,2) - val_b_RW1)/val_k_RW1;

% Calculating MOI of the RW
MOI_RW1 = (torque_const_RW * mean(Data_RW1(:,3)))/abs(p(1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Reaction Wheel Run with Torque = 10
[Time_RW2,Data_RW2] = LoadData_RWHEEL('Lab 3 Data/10-7-25_RWHEEL_T10.csv');
p = polyfit(Time_RW2,Data_RW2(:,2),1);
val_b_RW2 = p(2);% This is the bias in the RW
val_k_RW2 = p(1);% This is the scale factor for calibration
Polyline_RW2 = polyval(p,Time_RW2);

Calibrated_Data_RW2 = (Data_RW2(:,2) - val_b_RW2)/val_k_RW2;

% Calculating MOI of the RW
MOI_RW2 = (torque_const_RW * mean(Data_RW2(:,3)))/abs(p(1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Reaction Wheel Run with Torque = 15
[Time_RW3,Data_RW3] = LoadData_RWHEEL('Lab 3 Data/10-7-25_RWHEEL_T15.csv');
p = polyfit(Time_RW3,Data_RW3(:,2),1);
val_b_RW3 = p(2);% This is the bias in the RW
val_k_RW3 = p(1);% This is the scale factor for calibration
Polyline_RW3 = polyval(p,Time_RW3);

Calibrated_Data_RW3 = (Data_RW3(:,2) - val_b_RW3)/val_k_RW3;

% Calculating MOI of the RW
MOI_RW3 = (torque_const_RW * mean(Data_RW3(:,3)))/abs(p(1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loading Reaction Wheel Run with Torque = 20
[Time_RW4,Data_RW4] = LoadData_RWHEEL('Lab 3 Data/10-7-25_RWHEEL_T20.csv');
p = polyfit(Time_RW4,Data_RW4(:,2),1);
val_b_RW4 = p(2);% This is the bias in the RW
val_k_RW4 = p(1);% This is the scale factor for calibration
Polyline_RW4 = polyval(p,Time_RW4);

Calibrated_Data_RW4 = (Data_RW4(:,2) - val_b_RW4)/val_k_RW4;

% Calculating MOI of the RW
MOI_RW4 = (torque_const_RW * mean(Data_RW4(:,3)))/abs(p(1));


% Finding Standard Deviation and Mean MOI
[S, M] = std([MOI_RW1, MOI_RW2, MOI_RW3, MOI_RW4]); %mNm

% Calculating time to reach max angular velocity of the reaction wheel
maxAngularVel_RW = 4000; %rpm
maxAngularVel_RW_rad = 4000*(pi/180); %rad/s
applied_torque = 10e-4; %Nm
ang_accel_RW = applied_torque/(M/1000);
t_to_max_vel = maxAngularVel_RW_rad/ang_accel_RW;


%% BASE Plotting Section


%% CONTROL Plotting Section

if Control1Switch == 1
	figure(name='Control1 Kp=50 Kd=15')
	hold on
	grid on
	plot(Time_CONTROL1,wrapped_CONTROL1,LineWidth=1.8)
	plot(Time_CONTROL1,Data_CONTROL1(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=50, Kd=15 - CONTROL1',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL1(end)])
	print('Problem_3_3bv_CONTROL1','-dpng','-r600')
end

if Control2Switch == 1
	figure(name='Control2 Kp=100 Kd=15')
	hold on
	grid on
	plot(Time_CONTROL2,wrapped_CONTROL2,LineWidth=1.8)
	plot(Time_CONTROL2,Data_CONTROL2(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=100, Kd=15 - CONTROL2',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL2(end)])
	print('Problem_3_3bv_CONTROL2','-dpng','-r600')
end

if Control3Switch == 1
	figure(name='Control3 Kp=100 Kd=40')
	hold on
	grid on
	plot(Time_CONTROL3,wrapped_CONTROL3,LineWidth=1.8)
	plot(Time_CONTROL3,Data_CONTROL3(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=100, Kd=40 - CONTROL3',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL3(end)])
	print('Problem_3_3bv_CONTROL3','-dpng','-r600')
end

if Control4Switch == 1
	figure(name='Control4 Kp=100 Kd=-10')
	hold on
	grid on
	plot(Time_CONTROL4,wrapped_CONTROL4,LineWidth=1.8)
	plot(Time_CONTROL4,Data_CONTROL4(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=100, Kd=10 - CONTROL4',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL4(end)])
	print('Problem_3_3bv_CONTROL4','-dpng','-r600')
end

if Control5Switch == 1
	figure(name='Control5 Kp=150 Kd=15')
	hold on
	grid on
	plot(Time_CONTROL5,wrapped_CONTROL5,LineWidth=1.8)
	plot(Time_CONTROL5,Data_CONTROL5(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=150, Kd=15 - CONTROL5',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL5(end)])
	print('Problem_3_3bv_CONTROL5','-dpng','-r600')
end

if Control6Switch == 1
	figure(name='Control6 Kp=200 Kd=0')
	hold on
	grid on
	plot(Time_CONTROL6,wrapped_CONTROL6,LineWidth=1.8)
	plot(Time_CONTROL6,Data_CONTROL6(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=200, Kd=0 - CONTROL6',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL6(end)])
	print('Problem_3_3biii_CONTROL6','-dpng','-r600')
end

if ControlChosenSwitch == 1
	figure(name='Control Chosen Kp=3.48 Kd=16.6')
	hold on
	grid on
	plot(Time_CONTROL_Chosen, Data_CONTROL_Chosen(:,2)+0.06,LineWidth=1.8)
	plot(Time_CONTROL_Chosen, Data_CONTROL_Chosen(:,1),'r',LineWidth=1)
	title('S/C Control Angular Position w/ Kp=3.48, Kd=16.6 - Calibrated Control',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Position [rad]',FontSize=14)
	legend('S/C Angular Position','Reference Set Point',Location='northeast')
	hold off
	xlim([0,Time_CONTROL_Chosen(end)])
	print('Problem_3_3d_CONTROL_Chosen','-dpng','-r600')
end

%% GYRO Plotting Section

if Gyro1Switch == 1
	%figure(name='Gyro 1 Calibration')
	%axis square
	%hold on
	%grid on
	%scatter(Data_GYRO1(:,2),Data_GYRO1(:,1),'.');
	%plot([min_GYRO1,max_GYRO1],[maxPolyline_GYRO1,minPolyline_GYRO1],'--r',LineWidth=2)
	%yline(val_b_GYRO1,'--');
	%hold off

	% Plot for Problem 3.1 part c.4i
	figure(name='Calibrated GYRO1 Time Plot')
	hold on
	grid on
	plot(Time_GYRO1,Data_GYRO1(:,2),'--',LineWidth=1.8)
	plot(Time_GYRO1,Calibrated_Data_GYRO1,'r',LineWidth=0.8)
	title('Calibrated Angular Velocity vs Encoder Data - Gyro1',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data',Location='northwest')
	hold off
	xlim([0,Time_GYRO1(end)])
	print('Problem_3_1c4i_GYRO1','-dpng','-r600')

	% Plot for Problem 3.1 part c.4ii
	figure(name ='Error_GYRO1 in Angular Rate');
	hold on
	grid on
	plot(Time_GYRO1,angular_rate_error_GYRO1,LineWidth=0.8);
	xlabel('Time [s]',FontSize=14);
	ylabel('Angular Rate Error [rad/s]',FontSize=14);
	title('Time History of Angular Rate Error - Gyro1',FontSize=18);
	hold off
	xlim([0,Time_GYRO1(end)])
	print('Problem_3_1c4ii_GYRO1','-dpng','-r600')

	% Plot for Problem 3.1 part c.4iii
	figure(name='Gyro1_Position Measure Time Plot');
	hold on;
	grid on
	plot(Time_GYRO1, theta_GYRO1,'--',LineWidth=1.8);
	plot(Time_GYRO1, phi_GYRO1,'r',LineWidth=0.8);
	ylabel('Angular Position [rad]',FontSize=14);
	xlabel('Time [s]',FontSize=14);
	title('Time History of True Angular Position(Encoder) vs Measured Angular Position - Gyro1',FontSize=18)
	legend('True Angular Position','Measured Angular Position',Location='southeast');
	hold off
	xlim([0,Time_GYRO1(end)])
	print('Problem_3_1c4iii_GYRO1','-dpng','-r600')

	% Plot for Problem 3.1 part c.4iv
	figure(name ='Error_GYRO1 in Angular Rate');
	hold on
	grid on
	plot(Time_GYRO1,angular_position_error_GYRO1,LineWidth=0.8);
	xlabel('Time [s]',FontSize=14);
	ylabel('Angular Position Error [rad]',FontSize=14);
	title('Time History of Angular Position Error - Gyro1',FontSize=18);
	hold off
	xlim([0,Time_GYRO1(end)])
	print('Problem_3_1c4iv_GYRO1','-dpng','-r600')

	% Plot for Problem 3.1 part c.4v
	figure(name='Gyro1 Rate Error vs Position Error')
	axis square
	hold on
	grid on
	plot(Data_GYRO1(:,2),angular_position_error_GYRO1,LineWidth=0.8)
	xlabel('Encoder Angular Rate [rad/s]',FontSize=14);
	ylabel('Angular Position Error [rad]',FontSize=14);
	title('Angular Position Error as a Function of Encoder Angular Rate - Gyro1',FontSize=18);
	hold off
	print('Problem_3_1c4v_GYRO1','-dpng','-r600')
end

if Gyro2Switch == 1
	% Plot for Problem 3.1 part c.1.1
	figure(name='Raw Unflipped Gyro 2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,2),'--',LineWidth=2)
	plot(Time_GYRO2,Data_GYRO2(:,1),'r',LineWidth=0.8)
	title('Raw Gyro Data vs Encoder Data Before Calibration - Gyro2',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Encoder Data','Raw Gyro Data',Location='southeast')
	hold off
	xlim([0,Time_GYRO2(end)])
	print('Problem_3_1c1-1','-dpng','-r600')

	% Plot for Problem 3.1 part c.1.2
	figure(name='Raw Flipped Gyro 2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,2),'--',LineWidth=2)
	plot(Time_GYRO2,-Data_GYRO2(:,1),'r',LineWidth=0.8)
	title('Raw Gyro Data Flipped vs Encoder Data Before Calibration - Gyro2',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Encoder Data','Raw Gyro Data',Location='southeast')
	hold off
	xlim([0,Time_GYRO2(end)])
	print('Problem_3_1c1-2','-dpng','-r600')

	% Plot for Problem 3.1 part c.2
	figure(name='Gyro 2 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO2(:,2),Data_GYRO2(:,1),'.');
	plot([min_GYRO2,max_GYRO2],[maxPolyline_GYRO2,minPolyline_GYRO2],'--r',LineWidth=2)
	yline(val_b_GYRO2,'--');
	title('Gyro Calibration Data - Gyro2',FontSize=18)
	xlabel('Encoder Rate Measurement [rad/s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Data','Adjusted Scale Factor \itK','Bias \itb')
	hold off
	print('Problem_3_1c2-1','-dpng','-r600')
	
	% Plot for Problem 3.1 part c.2
	figure(name='Calibrated Gyro2 Time Plot')
	hold on
	grid on
	plot(Time_GYRO2,Data_GYRO2(:,2),'--',LineWidth=2)
	plot(Time_GYRO2,Calibrated_Data_GYRO2,'r',LineWidth=0.8)
	title('Calibrated Angular Velocity vs Encoder Data - Gyro2',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Outpu [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data',Location='southeast')
	hold off
	xlim([0,Time_GYRO2(end)])
	print('Problem_3_1c2-2','-dpng','-r600')
end

if Gyro3Switch == 1
	%figure(name='Gyro 3 Calibration')
	%axis square
	%hold on
	%grid on
	%scatter(Data_GYRO3(:,2),Data_GYRO3(:,1),'.');
	%plot([min_GYRO3,max_GYRO3],[maxPolyline_GYRO3,minPolyline_GYRO3],'--r',LineWidth=2)
	%yline(val_b_GYRO3,'--');
	%hold off
	
	% Plot for Problem 3.1 part c.4i
	figure(name='Calibrated GYRO3 Time Plot')
	hold on
	grid on
	plot(Time_GYRO3,Data_GYRO3(:,2),'--',LineWidth=1.8)
	plot(Time_GYRO3,Calibrated_Data_GYRO3,'r',LineWidth=0.8)
	title('Calibrated Angular Velocity vs Encoder Data - Gyro3',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Output [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data',Location='southeast')
	hold off
	xlim([0,Time_GYRO3(end)])
	print('Problem_3_1c4i_GYRO3','-dpng','-r600')

	% Plot for Problem 3.1 part c.4ii
	figure(name ='Error_GYRO3 in Angular Rate');
	hold on
	grid on
	plot(Time_GYRO3,angular_rate_error_GYRO3,LineWidth=0.8);
	xlabel('Time [s]',FontSize=14);
	ylabel('Angular Rate Error [rad/s]',FontSize=14);
	title('Time History of Angular Rate Error - Gyro3',FontSize=18);
	hold off
	xlim([0,Time_GYRO3(end)])
	print('Problem_3_1c4ii_GYRO3','-dpng','-r600')

	% Plot for Problem 3.1 part c.4iii
	figure(name='Gyro3_Position Measure Time Plot');
	grid on
	hold on;
	plot(Time_GYRO3, theta_GYRO3,'--',LineWidth=1.8);
	plot(Time_GYRO3, phi_GYRO3,'r',LineWidth=0.8);
	ylabel('Angular Position [rad]',FontSize=14);
	xlabel('Time [s]',FontSize=14);
	title('Time History of True Angular Position(Encoder) vs Measured Angular Position - Gyro3',FontSize=18)
	legend('True Angular Position','Measured Angular Position',Location='southeast');
	hold off
	xlim([0,Time_GYRO3(end)])
	print('Problem_3_1c4iii_GYRO3','-dpng','-r600')

	% Plot for Problem 3.1 part c.4iv
	figure(name ='Error_GYRO3 in Angular Rate');
	hold on
	grid on
	plot(Time_GYRO3,angular_position_error_GYRO3,LineWidth=0.8);
	xlabel('Time [s]',FontSize=14);
	ylabel('Angular Position Error [rad]',FontSize=14);
	title('Time History of Angular Position Error - Gyro3',FontSize=18);
	hold off
	xlim([0,Time_GYRO3(end)])
	print('Problem_3_1c4iv_GYRO3','-dpng','-r600')

	% Plot for Problem 3.1 part c.4v
	figure(name='Gyro3 Rate Error vs Position Error')
	axis square
	hold on
	grid on
	plot(Data_GYRO3(:,2),angular_position_error_GYRO3,LineWidth=0.8)
	xlabel('Encoder Angular Rate [rad/s]',FontSize=14);
	ylabel('Angular Position Error [rad]',FontSize=14);
	title('Angular Position Error as a Function of Encoder Angular Rate - Gyro3',FontSize=18);
	hold off
	print('Problem_3_1c4v_GYRO3','-dpng','-r600')
end

if Gyro4Switch == 1
	figure(name='Gyro 4 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO4(:,2),Data_GYRO4(:,1),'.');
	plot([min_GYRO4,max_GYRO4],[maxPolyline_GYRO4,minPolyline_GYRO4],'--r',LineWidth=2)
	yline(val_b_GYRO4,'--');
	hold off
	
	figure(name='Gyro 4 Time Plot')
	hold on
	grid on
	plot(Time_GYRO4,Data_GYRO4(:,2))
	plot(Time_GYRO4,Calibrated_Data_GYRO4)
	hold off
end

if Gyro5Switch == 1
	figure(name='Gyro 5 Calibration')
	axis square
	hold on
	grid on
	scatter(Data_GYRO5(:,2),Data_GYRO5(:,1),'.');
	plot([min_GYRO5,max_GYRO5],[maxPolyline_GYRO5,minPolyline_GYRO5],'--r',LineWidth=2)
	yline(val_b_GYRO5,'--');
	hold off
	
	figure(name='Gyro 5 Time Plot')
	hold on
	grid on
	plot(Time_GYRO5,Data_GYRO5(:,2))
	plot(Time_GYRO5,Calibrated_Data_GYRO5)
	hold off
end

if Gyro6Switch == 1
	%figure(name='Gyro 6 Calibration')
	%axis square
	%hold on
	%grid on
	%scatter(Data_GYRO6(:,2),Data_GYRO6(:,1),'.');
	%plot([min_GYRO6,max_GYRO6],[maxPolyline_GYRO6,minPolyline_GYRO6],'--r',LineWidth=2)
	%yline(val_b_GYRO6,'--');
	%hold off

	% Plot for Problem 3.1 part b.1
	figure(name='Gyro 6 Time Plot')
	hold on
	grid on
	plot(Time_GYRO6,Data_GYRO6(:,2),'--',LineWidth=1.8)
	plot(Time_GYRO6,Calibrated_Data_GYRO6,'r',LineWidth=0.8)
	title('Calibrated MANUAL Gyro Data vs Encoder Data - Gyro6',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data',Location='southwest')
	hold off
	xlim([0,Time_GYRO6(end)])
	print('Problem_3_1b-1','-dpng','-r600')

	% Plot for Problem 3.1 part b.2
	figure(name='Gyro 6 Time Plot')
	hold on
	grid on
	plot(Time_GYRO6,Data_GYRO6(:,2),'--',LineWidth=1.8)
	plot(Time_GYRO6,Data_GYRO6(:,1),'r',LineWidth=0.8)
	title('Raw MANUAL Gyro Data vs Encoder Data - Gyro6',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Raw Gyro Data',Location='southwest')
	hold off
	xlim([0,10]);
	print('Problem_3_1b-2','-dpng','-r600')

	% Plot for Problem 3.1 part b.3
	figure(name='Gyro 6 Time Plot')
	hold on
	grid on
	plot(Time_GYRO6,Data_GYRO6(:,2),'--',LineWidth=1.8)
	plot(Time_GYRO6,Calibrated_Data_GYRO6,'r',LineWidth=0.8)
	title('Calibrated MANUAL Gyro Data vs Encoder Data - Gyro6',FontSize=18)
	xlabel('Time [s]',FontSize=14)
	ylabel('Gyro Data [rad/s]',FontSize=14)
	legend('Encoder Data','Calibrated Gyro Data',Location='southwest')
	hold off
	xlim([0,10]);
	print('Problem_3_1b-3','-dpng','-r600')
end

%% REACTION WHEEL Plotting Section

if Rwheel1Switch == 1
    figure(name='Angular Velocity Reaction Wheel T5')
    hold on
    title('Angular Velocity Reaction Wheel T5',FontSize=18)
    scatter(Time_RW1,Data_RW1(:,2),'.');
    plot(Time_RW1, Polyline_RW1, 'r-', 'LineWidth', 1.8);
    legend('Raw Data', 'Line Fitted Data')
    xlabel('Time [s]') 
    ylabel('Angular Velocity [deg/s]')
    yline(val_b_RW1,'--');
    axis fill
    hold off
end

if Rwheel2Switch == 1
    figure(name='Angular Velocity Reaction Wheel T10')
    hold on 
    title('Angular Velocity Reaction Wheel T10',FontSize=18)
    scatter(Time_RW2,Data_RW2(:,2),'.');
    plot(Time_RW2, Polyline_RW2, 'r-', 'LineWidth', 1.8);
    legend('Raw Data', 'Line Fitted Data')
    xlabel('Time [s]')
    ylabel('Angular Velocity [deg/s]')
    yline(val_b_RW2,'--');
    axis fill
    hold off
end

if Rwheel3Switch == 1
    figure(name='Angular Velocity Reaction Wheel T15')
    hold on 
    title('Angular Velocity Reaction Wheel T15',FontSize=18)
    scatter(Time_RW3,Data_RW3(:,2),'.');
    plot(Time_RW3, Polyline_RW3, 'r-', 'LineWidth', 1.8);
    legend('Raw Data', 'Line Fitted Data')
    xlabel('Time [s]')
    ylabel('Angular Velocity [deg/s]')
    yline(val_b_RW3,'--');
    axis fill
    hold off
end

if Rwheel4Switch == 1
    figure(name='Angular Velocity Reaction Wheel T20')
    hold on 
    title('Angular Velocity Reaction Wheel T20',FontSize=18)
    scatter(Time_RW4,Data_RW4(:,2),'.');
    plot(Time_RW4, Polyline_RW4, 'r-', 'LineWidth', 1.8);
    legend('Raw Data', 'Line Fitted Data')
    xlabel('Time [s]')
    ylabel('Angular Velocity [deg/s]')
    yline(val_b_RW4,'--');
    axis tight
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
%data = unique(data, 'rows'); % Remove duplicate entries
Time = data(:, 1); % Extracting the time column
Data = data(:, 2:end); % Extracting the data columns

end

function [Time,Data] = LoadData_GYRO(fileName)
% LoadData is a function that brings in the .csv files

% Read in Data
data = readtable(fileName);
data = table2array(data);
data(:,3) = data(:,3)/60 * 2 * pi; % Convert ENCODER from RPM to rad/s
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