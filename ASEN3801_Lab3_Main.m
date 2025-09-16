% ASEN 2803 Rotary Arm Lab 3
% Authors: Brooklyn O'Laughlin, Andrew Gonzalez, Riley Martin, Ashton Unzicker
% Date: 4/8/2024
% Blah
% House Cleaning
close all; clc; clear

% Call to Constants function
const = getConst();

%% Gains set for Question 3.3

% Proportional gain applied to the hub angle measurement (-2 to 50)
K1_rigid = [10 20 5 10 10 10];
% Derivative gain applied to the hub angle rate (-1.5 to 1.5)
K3_rigid = [0 0 0 1 -1 -0.5];
B = 0;

omega = (K1_rigid * const.K_g * const.K_m) ./ (const.J * const.R_m);
zeta = (B * const.R_m + const.K_g ^ 2 * const.K_m ^ 2 + K3_rigid * const.K_g * const.K_m) ./ (2 * sqrt(K1_rigid * const.K_g * const.K_m * const.J * const.R_m));
time = linspace(0,10,100000);
U = 0.5 * ones(1,length(time));

figure(1)
hold on
for i = 1:length(K1_rigid)

	% Call to armEOM function
	[xs,t] = armEOM(omega(i),zeta(i),U,time);
	
	plot(t,xs,LineWidth=1.5)
	title('Path to 0.5 Radians Through Different KP and KD Values');
	ylabel('Theta (rad)');
	xlabel('Time (s)');
end

% Plotting for Question 3.3
figure(1)
yline(0.5)
ax = gca;
ax.FontSize = 14;
legend( "KP = " + K1_rigid(1) + ", KD = " + K3_rigid(1), ...
		"KP = " + K1_rigid(2) + ", KD = " + K3_rigid(2), ...
		"KP = " + K1_rigid(3) + ", KD = " + K3_rigid(3), ...
		"KP = " + K1_rigid(4) + ", KD = " + K3_rigid(4), ...
		"KP = " + K1_rigid(5) + ", KD = " + K3_rigid(5), ...
		"KP = " + K1_rigid(6) + ", KD = " + K3_rigid(6))

%% Gains set for Question 4.3 b)

% Proportional gain applied to the hub angle measurement (-2 to 50)
K1_rigid = 9;
% Derivative gain applied to the hub angle rate (-1.5 to 1.5)
K3_rigid = 0.5;
B = 0;

omega = (K1_rigid * const.K_g * const.K_m) ./ (const.J * const.R_m);
zeta = (B * const.R_m + const.K_g ^ 2 * const.K_m ^ 2 + K3_rigid * const.K_g * const.K_m) ./ (2 * sqrt(K1_rigid * const.K_g * const.K_m * const.J * const.R_m));
time = linspace(0,10,100000);
U = [0.5 * ones(1,length(time)/2),...
	-0.5 * ones(1,length(time)/2)];

% Call to armEOM function
[xs,t] = armEOM(omega,zeta,U,time);

% Plotting for Question 4.3
figure(2)
hold on
plot(t,xs,LineWidth=1.5)
plot(time,U,'m',LineWidth=1.5)
yline(-0.525,'--b',LineWidth=1.5)
yline(-0.475,'--b',LineWidth=1.5)
yline(-0.6,'--r',LineWidth=1.5)
yline(-0.4,'--r',LineWidth=1.5)
yline(0)
%xline(6)
title('Control Gains Set to Meet Requirements');
ax = gca;
ax.FontSize = 14;
ylabel('Theta (rad)');
xlabel('Time (s)');
legend("KP = " + K1_rigid + ", KD = " + K3_rigid, 'Set Point', '5% Deviation','','20% Deviation')
ylim([-0.65,0.65])

%% LabView Plot Overlay Section 4.3 c)

time = linspace(0,20,100000);
U = [0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4),...
	0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4)];

% Cleaning LabView Data
LabView = load("Group_05");
reset_row = find(LabView(:,6)>0,1);
LabView(1:reset_row,:) = [];
LabView(:,1) = (LabView(:,1) - LabView(1,1)) / 1000;

figure(3)
hold on
plot(LabView(:,1),LabView(:,2),'r',LineWidth=2);
plot(time,U,'m',LineWidth=1.5)
yline(-0.525,'--b',LineWidth=1.5)
yline(-0.475,'--b',LineWidth=1.5)
yline(-0.6,'--k',LineWidth=1.5)
yline(-0.4,'--k',LineWidth=1.5)
yline(0.525,'--b',LineWidth=1.5)
yline(0.475,'--b',LineWidth=1.5)
yline(0.6,'--k',LineWidth=1.5)
yline(0.4,'--k',LineWidth=1.5)

yline(0)
title('LabView Data Plot (Matching Gains)');
ax = gca;
ax.FontSize = 14;
ylabel('Theta (rad)');
xlabel('Time (s)');
legend('LabView Data','Set Point','5% Deviation','','20% Deviation',Location='north')
xlim([2.5,12.5])
ylim([-0.65,0.65])

%% LabView Plot Overlay Section 4.3 d)

% Proportional gain applied to the hub angle measurement (-2 to 50)
K1_rigid = 9;
% Derivative gain applied to the hub angle rate (-1.5 to 1.5)
K3_rigid = 0.5; % 1.25 Corrected
B = 0;

omega = (K1_rigid * const.K_g * const.K_m) ./ (const.J * const.R_m);
zeta = (B * const.R_m + const.K_g ^ 2 * const.K_m ^ 2 + K3_rigid * const.K_g * const.K_m) ./ (2 * sqrt(K1_rigid * const.K_g * const.K_m * const.J * const.R_m));
time = linspace(0,20,100000);
U = [0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4),...
	0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4)];

% Cleaning LabView Data
LabView = load("Group_05");
reset_row = find(LabView(:,6)>0,1);
LabView(1:reset_row,:) = [];
LabView(:,1) = (LabView(:,1) - LabView(1,1)) / 1000;

% Call to armEOM function
[xs,t] = armEOM(omega,zeta,U,time);

figure(4)
hold on
plot(t,xs,LineWidth=2)
plot(LabView(:,1),LabView(:,2),'r',LineWidth=2);
plot(time,U,'m',LineWidth=1.5)
yline(-0.525,'--b',LineWidth=1.5)
yline(-0.475,'--b',LineWidth=1.5)
yline(-0.6,'--k',LineWidth=1.5)
yline(-0.4,'--k',LineWidth=1.5)
yline(0.525,'--b',LineWidth=1.5)
yline(0.475,'--b',LineWidth=1.5)
yline(0.6,'--k',LineWidth=1.5)
yline(0.4,'--k',LineWidth=1.5)

yline(0)
title('Overlay - Simulation vs LabView Data Plot');
ax = gca;
ax.FontSize = 14;
ylabel('Theta (rad)');
xlabel('Time (s)');
legend("KP = " + K1_rigid + ", KD = " + K3_rigid,'LabView Data','Set Point','5% Deviation','','20% Deviation',Location='north')
xlim([2.5,12.5])
ylim([-0.65,0.65])

%% Gains set for Question 5.2 b) Natural Dampening

% Proportional gain applied to the hub angle measurement (-2 to 50)
K1_rigid = 9;
% Derivative gain applied to the hub angle rate (-1.5 to 1.5)
K3_rigid = 0.5;
K3_rigid_DAMPED = 1.25;
B = 0;

omega = (K1_rigid * const.K_g * const.K_m) ./ (const.J * const.R_m);
zeta = (B * const.R_m + const.K_g ^ 2 * const.K_m ^ 2 + K3_rigid * const.K_g * const.K_m) ./ (2 * sqrt(K1_rigid * const.K_g * const.K_m * const.J * const.R_m));
omega_DAMPED = (K1_rigid * const.K_g * const.K_m) ./ (const.J * const.R_m);
zeta_DAMPED = (B * const.R_m + const.K_g ^ 2 * const.K_m ^ 2 + K3_rigid_DAMPED * const.K_g * const.K_m) ./ (2 * sqrt(K1_rigid * const.K_g * const.K_m * const.J * const.R_m));
time = linspace(0,20,100000);
U = [0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4),...
	0.5 * ones(1,length(time)/4),...
	-0.5 * ones(1,length(time)/4)];

% Cleaning LabView Data
LabView = load("Group_05");
reset_row = find(LabView(:,6)>0,1);
LabView(1:reset_row,:) = [];
LabView(:,1) = (LabView(:,1) - LabView(1,1)) / 1000;

% Call to armEOM function
[xs,t] = armEOM(omega,zeta,U,time);
[xs_DAMPED,t_DAMPED] = armEOM(omega_DAMPED,zeta_DAMPED,U,time);

figure(5)
hold on
plot(t,xs,LineWidth=2)
plot(t_DAMPED,xs_DAMPED,'g',LineWidth=3)
plot(LabView(:,1),LabView(:,2),'r',LineWidth=2);
plot(time,U,'m',LineWidth=1.5)
yline(-0.525,'--b',LineWidth=1.5)
yline(-0.475,'--b',LineWidth=1.5)
yline(-0.6,'--k',LineWidth=1.5)
yline(-0.4,'--k',LineWidth=1.5)
yline(0.525,'--b',LineWidth=1.5)
yline(0.475,'--b',LineWidth=1.5)
yline(0.6,'--k',LineWidth=1.5)
yline(0.4,'--k',LineWidth=1.5)

yline(0)
title('Overlay - Simulation vs Damped Simulation vs LabView Data Plot');
ax = gca;
ax.FontSize = 14;
ylabel('Theta (rad)');
xlabel('Time (s)');
legend("KP = " + K1_rigid + ", KD = " + K3_rigid,"KP = " + K1_rigid + ", KD = " + K3_rigid_DAMPED,'LabView Data','Set Point','5% Deviation','','20% Deviation',Location='north')
xlim([2.5,12.5])
ylim([-0.65,0.65])

%% Functions Section

function [xs,t] = armEOM(omega,zeta,U,time)
	% armEOM

	val = omega;
	d2 = 1;
	d1 = 2 * zeta * sqrt(omega);
	d0 = omega;
	den = [d2 d1 d0];
	sysTF = tf(val,den);
	
	% lsim Specific u(t) function
	[xs,t] = lsim(sysTF,U,time);
	
	% Step function
	%[Xs,tS] = step(sysTF);
	% figure()
	% plot(tS,Xs)
end

function const = getConst()
	% getConst

	% Base
	const.K_g = 33.3; % Gear Ratio [No Units]
	const.K_m = 0.0401; % Proportional motor constant [V/rad/s]
	const.R_m = 19.2; % Output motor resistance [ohms]
	
	% Rigid Arm
	const.J_hub = 0.0005; % Inertia hub [kg*m^2]
	const.J_load = 0.0015; % Inertia due to added load [kg*m^2]
	const.J_extra = (0.2 * 0.2794 ^ 2); % Inertia of extra mass [kg*m^2]
	const.J = const.J_hub + const.J_load + const.J_extra; % Total intertia [kg*m^2]
	
	% Flexible Link
	const.J_hub1 = 0.0005; % [kg*m^2]
	const.L = 0.45; % Link Length [m]
	const.M_arm = 0.06; %kg Link mass of ruler [kg]
	const.J_arm = (const.M_arm * const.L ^ 2) / 3; % [kg*m^2]
	const.M_tip = 0.05; % [kg]
	const.J_M = const.M_tip * const.L ^ 2; % [kg*m^2]
	const.fc = 1.8; % Natural frequency [Hz]
	const.J_L = const.J_arm + const.J_M; % Intertia of Flex Link + Mass
	const.Karm = (2 * pi * const.fc)^2 * const.J_L; %Flexible Link Stiffness (Kstiff)

	%{
	Parking for Flexible Arm Gains (unused at the moment)
	% Flexible Arm Gains
	% Proportional gain applied to the hub angle measurement (0 to 20)
	K1_flex = [0 0 0 0 0 0];
	% Proportional gain applied to the tip sensor measurement (-50 to 0)
	K2_flex = [0 0 0 0 0 0];
	% Derivative gain applied to the hub angle rate (0 to 1.5)
	K3_flex = [0 0 0 0 0 0];
	% Derivative gain applied to the tip sensor rate (0.1 to 0.3)
	K4_flex = [0 0 0 0 0 0];
	%}

end