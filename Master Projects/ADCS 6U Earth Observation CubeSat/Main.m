%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  SPACECRAFT ATTITUDE DYNAMICS                           %
%                    Academic year 2023/2024                              %
%                    M.Sc. Space Engineering                              %
%                     Politecnico di Milano                               %
%                                                                         %
%         Simulation of a Spacecraft's Attitude Dynamics and Control      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Group N.26
% Project N.217 
% Platform: Cubesat 6U
% Attitude Parameters: Quaternions
% Mandatory Sensor: Magnetic Field Sensor
% Actuators: 3 Reaction wheels
%
% Authors: - Casiero Alessia
%          - Fiume Elisa
%          - Gallo Emanuele
%          - Marotta Arianna

%%
close all
clear
clc

% Path to the functions:
addpath(strcat(pwd,'/Functions'));
addpath(strcat(pwd,'/Plots'));
%% Mass Properties 
% Developed by designing a model of the Cubesat 6U assigned on SolidWorks

% Inertia matrix computed for the Cubesat with 4 deployed solar panels [kg*m^2]
Ix = 0.13;
Iy = 0.16;
Iz = 0.18;
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % [kg*m^2]
I_inv = inv(I);
r_CG = [0.02 0 0]; % Baricenter position

%% Orbit
% A Sun-synchronus orbit has been deployed 

R_E = astroConstants(23)*1000; % Earth Radius [m]
G = 6.67*10^-11; % Gravitational constant
e = 0.00149; % Eccentricity
i = deg2rad(97.525); % Inclination [rad]
OM = deg2rad(6.934*180/12); % Right Ascension of Ascending Node [rad]
om = deg2rad(287.373); % Argument of Perigee [rad]
th = deg2rad(72.588); % True Anomaly [rad]
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
a = (509.4*1000+R_E+529.9*1000+R_E)/2; % Semi-major axis [m]
w_E = 7.291597763887421e-05; % Earth Angular Velocity [rad/s]
Mt = 5.97219e24; 
n = sqrt(G*Mt/a^3); 

% Sun Properties 
T_earth = 365*24*3600; % Earth's Orbital period [s]
n_sun = 2*pi/T_earth; % Sun's mean Angular Velocity [rad]
epsilon = deg2rad(23.45); % Ecliptic [rad]
R_sun = astroConstants(2)*1e3; % Distance Sun - Earth % [m]

%% Dynamics
% Initial conditions:
w_x0 = deg2rad(-5); % [rad/s]
w_y0 = deg2rad(3); % [rad/s]
w_z0 = deg2rad(+2); % [rad/s]
w0 = [w_x0 w_y0 w_z0];

theta0 = 0;
w_LN = [0 0 n];

A_LN0 = [cos(n*0) sin(n*0) 0; -sin(n*0) cos(n*0) 0; 0 0 1];
A_BN0 = A_LN0;

j_B = [0.01, 0.05, 0.01]';

%% SENSORS
% Contains the mandatory sensor (Magnetic Field sensor) and two additional
% ones to improve the performances (Sun Sensor & Gyroscopes)

% Magnetic Field Sensor
sensors.mag.lin = 0.1; % [%FS] - Linearity Error
sensors.mag.rep = 0.1; % [%FS] - Repeatibility Error
sensors.mag.hys = 0.06; % [%FS] - Hysteresis Error
sensors.mag.FS = 6; % [Gauss] - Full-Scale
sensors.mag.res = 120e-6; % [Gauss] - Resolution
sensors.mag.bandwidth = 5e+6; % [Hz]
sensors.mag.TS = 1/sensors.mag.bandwidth;
sensors.mag.mys1 = deg2rad(0.01); % Non orthogonality wrt principal axis of measurement
sensors.mag.mys2 = deg2rad(1); % Non orthogonality wrt secundary axes of measurement
sensors.mag.accuracy = 15; 

% Gyroscopes
sensors.gyro.Fs = 262; % [Hz] - Sampling frequency
sensors.gyro.Ts = 1/sensors.gyro.Fs; % [s] - Sampling time
sensors.gyro.ARW = deg2rad(0.15)/sqrt(3600); % [rad/(s^(1/2))]
sensors.gyro.RRW = deg2rad(0.0003)/(3600)^3/2; % [rad/s^3/2]
sensors.gyro.sigma_n = sensors.gyro.ARW/sqrt(sensors.gyro.Ts); 
sensors.gyro.sigma_b = sensors.gyro.RRW/sqrt(sensors.gyro.Ts);
sensors.gyro.RRW0 = sensors.gyro.sigma_b^2*sensors.gyro.Ts;
sensors.gyro.mys = 1e-3; % [rad]
sensors.gyro.FS = deg2rad(400); % [rad/s]
sensors.gyro.resolution = deg2rad(0.22)/3600; % [rad/s]

% Sun Sensor
sensors.SS.accuracy = 0.5; % Accuracy of sensor
FOV = deg2rad(114); % Field of View

% ROTATION MATRICES & Orientation of the six Sun Sensors (normal wrt z-axis
% of Body Frame)
theta = [0 pi -pi/2 pi/2 -pi/2 pi/2];
A_ssb = zeros(3,3,6);

for j = 1:6
    if j<4
        A_ssb(:,:,j) = [ 1       0          0; ...                          % Rotation Matrix Body to Sensor i
                         0  cos(theta(j))  sin(theta(j)); ...
                         0  -sin(theta(j)) cos(theta(j))];
    else 
        A_ssb(:,:,j) = [cos(theta(j))           0       sin(theta(j)); ...  % Rotation Matrix Body to Sensor i
                            0                   1              0; ...
                       -sin(theta(j))           0      cos(theta(j))];
    end
end

%% ACTUATORS
% Contains the mandatory actuators (3 Reaction Wheels) and the additional
% Magnetorquer that works in case of saturation of the RWs

% Reaction Wheels -TA6494 Series
actuators.RW.Msat = 0.02; % [Nms]
actuators.RW.hsat = 0.3; % [Nms]
actuators.RW.wsat = 5000/60*2*pi; % [rad/s]
actuators.RW.Ir_singleRW = 7e-4; % [kg*m^2]
actuators.RW.MountingError = 0.1; % [% of inertia moment of the single actuator]
actuators.RW.A = eye(3,3); % Configuration matrix of actuators
actuators.RW.Ir = actuators.RW.Ir_singleRW*actuators.RW.A + ...
                  eye(3,3)*actuators.RW.MountingError/100*actuators.RW.Ir_singleRW; 
                  % Inertia matrix of the RWs' configuration (considering mounting error)

%% ENVIRONMENT - DISTURBING TORQUES
% Earth's Magnetic Field
environment.R = R_E; 
settings.N=13;
[WMM, environment.WMM.K, environment.WMM.g, environment.WMM.h] = IGRF_coeffs(settings.N);

% Earth's Magnetic Field in developed on Simulink due to a Matlab function
% exploiting the IGRF up to the 13th order, since we're orbiting on a LEO

% Solar Radiation Pressure
F_e = 1358+600+150; % Solar Radiation Intensity [W/m^2]
c = astroConstants(5)*1e3; % Speed of light [m/s]
P = F_e/c; % Solar Pressure [Pa]

% Geometry
% Normal vectors: Main S/C
n1 = [-1,0,0];
n2 = [0,1,0];
n3 = [0,0,1];
n4 = -n1;
n5 = -n2;
n6 = -n3;
% Normal vectors: Solar Arrays
n7 = n1;
n8 = n1;
n9 = n1;
n10 = n1;
n11 = n4;
n12 = n4;
n13 = n4;
n14 = n4;
N_B = [n1; n2; n3; n4; n5; n6; n7; n8; n9; n10;n11;n12;n13;n14];

% Areas: Main S/C [m^2]
A1 = 2*10^-2;
A2 = 3*10^-2;
A3 = 6*10^-2;
A4 = A1;
A5 = A2;
A6 = A3;
% Areas: Solar Arrays [m^2]
A7 = A3;
A8 = A7;
A9 = A2;
A10 = A2;
A11 = A7;
A12 = A7;
A13 = A9;
A14 = A9;

A = [A1; A2; A3; A4; A5; A6; A7; A8; A9; A10; A11; A12; A13; A14];
 
% Distance from the center of mass: Main S/C [m]
U = 0.1; % [m]
rF1 = [-1.5,0,0]*U;
rF2 = [0,1,0]*U;
rF3 = [0,0,0.5]*U;
rF4 = [1.5,0,0]*U;
rF5 = [0,-1,0]*U;
rF6 = [0,0,-0.5]*U;
% Distance from the center of mass: Solar Arrays [m]
rF7 = [1.5,0,2]*U;
rF8 = [1.5,0,-2]*U;
rF9 = [1.5,2.5,0]*U;
rF10 = [1.5,-2.5,0]*U;
rF11 = rF7;
rF12 = rF8;
rF13 = rF9;
rF14 = rF10;

rF = [rF1; rF2; rF3; rF4; rF5; rF6; rF7; rF8; rF9; rF10;rF11;rF12;rF13;rF14];

q0 = [0 0 0 1];

%% CONTROL
% State Observer
sensors.Observer.Lw = 0.05; % L matrix for w observation
sensors.Observer.Lq = 0.6; % L matrix for attitude observation

% Detumbling 
control.detumbling.MaxTime = 1000; % [s]
control.detumbling.wd = [0, 0, 0]; % [rad/s]
control.detumbling.kd = 0.03;

% Slew Manoeuvre (Non-linear Control)
control.slewManeuvre.MaxTime = 4000; % [s]
control.slewManeuvre.wd = [0, 0, 0]; % [rad/s]
control.slewManeuver.k1 = 0.4;
control.slewManeuver.k2 = 0.06;

% Trajectory Tracking  
control.tracking.wd = [0, 0, n]; % [rad/s]
control.tracking.k1 = 0.5;
control.tracking.k2 = 0.07;

%% Simulation
Tsample = 0.1; % Chosen time sample to simulate the sensors discrete sampling 
out = sim('Project26.slx', 'FixedStep', '0.01', 'Solver', 'ode4', 'StopTime', '4*pi/n');

%% Plots
choice = input("Do you want to visualize the plots? Insert 1 (True) or 0 (False)\n");
if choice
    run plots.m
end
