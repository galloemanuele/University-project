% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2, Exercise 3
% Author: Emanuele Gallo

% Prepare the workspace:
clc; clearvars; cspice_kclear; close all

% Load kernels:
cspice_furnsh('assignment02.tm');

% Allocate the initial state of the lunar orbiter in Moon Centred Inertial 
% (MCI) reference system:
rr0_orb = [4307.844185282820, -1317.980749248651, 2109.210101634011]';  % Initial position
vv0_orb = [-0.110997301537882, -0.509392750828585, 0.815198807994189]'; % Initial velocity
S0_orb = [rr0_orb; vv0_orb];                                            % Initial state

% Initial time:
et0_str =  '2024-11-18T16:30:00.000';  % String format
et0 = cspice_str2et(et0_str);          % Ephemeris Time (ET) format

% Final time:
etf_str = '2024-11-18T20:30:00.000';    % String format
etf = cspice_str2et(etf_str);           % Ephemeris Time (ET) format

% Measurement noise:
sig_p = 0.1;   % [km]: Measurememt noise

% Initial covariance matrix:
P0 = diag([10, 1, 1, 0.001, 0.001, 0.001, 0.00001, 0.00001]); % [km^2, km^2/s^2, rad^2]

% Moon useful constants:
par.M_Radii = cspice_bodvrd('Moon', 'RADII', 3);  % [km]: Moon radii
par.mu_M = cspice_bodvrd('Moon', 'GM', 1);        % [km^3/s^2]: Moon gravitational constant

% Initial guess coordinates:
lat0 = 78*cspice_rpd;                 % [rad]: Latitude
lon0 = 15*cspice_rpd;                 % [rad]: Longitude
h0 = 0;                               % [km]: Altitude
rho0 = par.M_Radii(1) + h0;           % [km]: Initial range
par.LL0_land = [rho0; lat0; lon0];    % [km, rad, rad]: Initial latitudinal coordinates

% Mask angle:
par.mask_angle = 0;   % [rad]

%% 3.1

% Build the time span vector:
time_step = 30;                             % [s]: time-step length
npoints = round((etf-et0)/time_step)+1;     % [-]: number of points in the interval
etspan = linspace(et0, etf, npoints);       % [s]: Vector of times

% Check the visibility windows and the relatve values of latitudinal
% coordinates:
[El_visible, vw, xx_orb_MCI] = visibility_windows(S0_orb, etspan, par);

% Check over the visibility over the entire time interval:
fprintf('-----------------------------------Visibility window---------------------------------------\n')
if isequal(length(vw), length(etspan))
    fprintf('The lander and the orbiter are in relative visibility for the entire time interval. \n\n')
else 
    fprintf('The lander and the orbiter are not in relative visibility for the entire time interval. \n\n')   
end

% Plot the evolution in time of the elevation:
figure;
hold on;
grid on;
time_hours = (vw - vw(1)) / 3600;
plot(time_hours, El_visible * cspice_dpr(), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Elevation');
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Elevation [deg]', 'Interpreter', 'latex', 'FontSize', 16);
title('Elevation profile over the visibility window (in MOONLANDER\_TOPO)', 'Interpreter', 'latex', 'FontSize', 11);
xline(time_hours(end), 'LineStyle', '--', 'HandleVisibility', 'off', 'Color', [0 1 1], 'LineWidth', 2);
yline(par.mask_angle * cspice_dpr(), 'r--', 'LineWidth', 2, 'HandleVisibility','off');
text(mean(time_hours), par.mask_angle * cspice_dpr() + 5, 'Minimum elevation', ...
     'Color', 'r', 'FontSize', 16, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
text(time_hours(end) + 0.3, +8, '$t_f$', 'Color', [0 1 1], ...
     'FontSize', 16, 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
ylim([-5, 90]);
grid on;
hold off;


%% 3.2

% Time evolution of the orbiter position vector:
xx_orb_MCI_exact = xx_orb_MCI;   % This was already computed within the visibility window function

% Time evolution of the relative range between the orbiter and the lander:
landerName = 'MOONLANDER';                                                      % Set the lander name
lander_pos_MCI = cspice_spkpos(landerName, etspan, 'J2000', 'NONE', 'MOON');    % Extract its position from kernels in J2000 coordinates with respect to Moon
rel_pos_MCI = xx_orb_MCI_exact(:, 1:3)' - lander_pos_MCI;                       % Calculate the relative position
rel_range = vecnorm(rel_pos_MCI);                                               % Calculate the range

% Exact solution for latitude and longitude:
lander_pos_IAUMOON = cspice_spkezr('MOONLANDER', etspan(1), 'IAU_MOON', 'NONE', 'MOON');  % Extract from kernels the lander position in IAU moon
[~, lon_real, lat_real] = cspice_reclat(lander_pos_IAUMOON(1:3));                         % Convert the lander position in latitudinal coordinates

% Simulate the measurements:
par.noise_cov = sig_p.^2 * eye(4);                                      % [km^2]: Noise covariance matrix
rng('swb2712')                                                          % Set the default algorithm to generate the random numbers
meas = mvnrnd([xx_orb_MCI_exact(:, 1:3), rel_range'], par.noise_cov);   % Run the mvnrnd function to obtain the measurements

% Calculate the errors:
err_pos = meas(:, 1:3) - xx_orb_MCI_exact(:, 1:3);  % Error on position [km]
err_range = meas(:, 4) - rel_range';                % Error on range [km/s]

% Plot the errors:
figure;
hold on;
grid on;
time_hours = (vw - vw(1))/3600;                                                         % Extract time in hours
plot(time_hours, err_pos(:, 1), 'LineWidth', 1, 'DisplayName', 'Position Error (X)');   % Plot error on first component of position
plot(time_hours, err_pos(:, 2), 'LineWidth', 1, 'DisplayName', 'Position Error (Y)');   % Plot error on second component of position
plot(time_hours, err_pos(:, 3), 'LineWidth', 1, 'DisplayName', 'Position Error (Z)');   % Plot error on third component of position
plot(time_hours, err_range, 'LineWidth', 1, 'DisplayName', 'Range Error');              % Plot range error
yline(3*sig_p, '--k', 'LineWidth', 1.2, 'DisplayName', '$+3 \sigma$ Threshold', 'HandleVisibility','off');        % Add +3 sigma lines for measurement noise
yline(-3*sig_p, '--k', 'LineWidth', 1.2, 'DisplayName', '$-3 \sigma$ Threshold', 'DisplayName', '$\pm 3\sigma$'); % Add -3 sigma lines for measurement noise
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Measurement Error [km]', 'Interpreter', 'latex', 'FontSize', 14);
title('Time evolution of the measurement errors @Moon J2000 frame', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
hold off;


%% 3.3

% Unscented Transform parameters:
par.alpha = 0.01;
par.beta = 2;

% Retrieve the initial state:
S0_orb_per = mvnrnd(S0_orb, P0(1:6, 1:6));  % Randomly perturb the initial condition

% Run the Unscented Kalman Filter (UKF):
[E_UKF_3, P_UKF_3] = UKF(S0_orb_per, P0(1:6, 1:6), etspan, meas(:, 1:3), par, '1'); % Only processing the position vector measurements

% Calculate errors:
err_pos = vecnorm(E_UKF_3(1:3, :) - xx_orb_MCI_exact(:, 1:3)');    % Norm of the position errors [km] (Case 1)
err_vel = vecnorm(E_UKF_3(4:6, :) - xx_orb_MCI_exact(:, 4:6)');    % Norm of the velocity errors [km/s] (Case 1)

% Calculate estimated covariance bounds (3-sigma):
est_cov_pos = 3 * sqrt(arrayfun(@(k) trace(P_UKF_3(1:3, 1:3, k)), 1:size(P_UKF_3, 3)));     % Covariance bound for the position [km] (Case 1)
est_cov_vel = 3 * sqrt(arrayfun(@(k) trace(P_UKF_3(4:6, 4:6, k)), 1:size(P_UKF_3, 3)));     % Covariance bound for the velocity [km/s] (Case 1)


% Plot the error on the position:
figure
semilogy(time_hours, err_pos, 'LineWidth', 1.5, 'Color',[0 0.4470 0.7410], 'DisplayName','$||\mathbf{r}_{\mathrm{err}}||$'); hold on
semilogy(time_hours, est_cov_pos, '--r', 'LineWidth', 1.5, 'DisplayName','$3\sigma_{pc}$');
hold off
grid on
xlabel('Time from $t_0$[h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Position Error [km]', 'Interpreter', 'latex', 'FontSize', 12);
title('Position Error and 3-Sigma Bound (time evolution in Moon J2000)', 'Interpreter', 'latex', 'FontSize', 10);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'Best');

% Plot the error on the velocity:
figure
semilogy(time_hours, err_vel, 'b', 'LineWidth', 1.5, 'Color',[0 0.4470 0.7410], 'DisplayName', '$||\mathbf{v}_{\mathrm{err}}||$'); hold on
semilogy(time_hours, est_cov_vel, '--r', 'LineWidth', 1.5,'DisplayName', '$3\sigma_{vc}$');
hold off 
grid on
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Velocity Error [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Velocity Error and 3-Sigma Bound (time evolution in Moon J2000)', 'Interpreter', 'latex', 'FontSize', 10);
legend('show', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'Best');

fprintf('3.3)---------------------------Estimate of SC state-------------------------------\n')
fprintf('-----------------------using its position vector measurements---------------------\n')

% Display lunar lander absolute state at initial time
fprintf('The estimated state at initial time is:\n');
disp(E_UKF_3(:, 1))
fprintf('\n')

% Display lunar lander absolute state at final time
fprintf('The estimated state at final time is:\n');
disp(E_UKF_3(:, end))
fprintf('With related covariance matrix:\n');
for i = 1:size(P_UKF_3, 1)
    fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P_UKF_3(i, :, end));
end
fprintf('\n\n\n')



%% 3.4

% Create the initial state as prescribed:
S0_ext_per = mvnrnd([S0_orb; par.LL0_land(2); par.LL0_land(3)], P0); % Perturbed initial state [km, km/s, rad]

% Run the Unscented Kalman Filter:
[E_UKF_4, P_UKF_4] = UKF(S0_ext_per, P0, etspan, meas, par, '2');

% Calculate errors:
err_pos4 = vecnorm(E_UKF_4(1:3, :) - xx_orb_MCI_exact(:, 1:3)');                                % Norm of the position errors [km]
err_vel4 = vecnorm(E_UKF_4(4:6, :) - xx_orb_MCI_exact(:, 4:6)');                                % Norm of the velocity errors [km/s]
err_lat = abs(angdiff(E_UKF_4(7, :), repmat(lat_real, 1, size(E_UKF_4, 2))));                   % Norm of the latitude errors [rad]
err_lon = abs(angdiff(E_UKF_4(8, :), repmat(lon_real, 1, size(E_UKF_4, 2))));                   % Norm of the longitude errors [rad]

% Calculate estimated covariance bounds (3-sigma):
est_cov_pos4 = 3 * sqrt(arrayfun(@(k) trace(P_UKF_4(1:3, 1:3, k)), 1:size(P_UKF_4, 3)));   % Covariance bound for the position [km]
est_cov_vel4 = 3 * sqrt(arrayfun(@(k) trace(P_UKF_4(4:6, 4:6, k)), 1:size(P_UKF_4, 3)));   % Covariance bound for the velocity [km/s]
est_cov_lon = 3 * sqrt(arrayfun(@(k) P_UKF_4(8, 8, k), 1:size(P_UKF_4, 3)));               % Covariance bound for the azimuth [rad]
est_cov_lat = 3 * sqrt(arrayfun(@(k) P_UKF_4(7, 7, k), 1:size(P_UKF_4, 3)));               % Covariance bound for the elevation [rad]

% Plot the error on the position:
figure
semilogy(time_hours, err_pos4, 'b', 'LineWidth', 1.5,'Color',[0 0.4470 0.7410]); hold on
semilogy(time_hours, est_cov_pos4, '--r', 'LineWidth', 1.5);
hold off
grid on
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Position Error [km]', 'Interpreter', 'latex', 'FontSize', 12);
title('Position Error and 3-Sigma Bound (time evolution in Moon J2000)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'$||\mathbf{r}_{\mathrm{err}}||$', '$3\sigma$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');

% Plot the error on the velocity:
figure
semilogy(time_hours, err_vel4, 'b', 'LineWidth', 1.5,'Color',[0 0.4470 0.7410]); hold on
semilogy(time_hours, est_cov_vel4, '--r', 'LineWidth', 1.5);
hold off
grid on
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Velocity Error [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Velocity Error and 3-Sigma Bound (time evolution in Moon J2000)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'$||\mathbf{v}_{\mathrm{err}}||$', '$3\sigma$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');

% Plot the error on the latitude:
figure
semilogy(time_hours, err_lat, 'b', 'LineWidth', 1.5,'Color',[0 0.4470 0.7410]); hold on
semilogy(time_hours, est_cov_lat, '--r', 'LineWidth', 1.5);
hold off
grid on
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Latitude Error [rad]', 'Interpreter', 'latex', 'FontSize', 12);
title('Latitude Error and 3-Sigma Bound (time evolution in MOONLANDER\_TOPO)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'Latitude error', '$3\sigma_{\mathrm{lat}}$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');

% Plot the error on the longitude:
figure
semilogy(time_hours, err_lon, 'b', 'LineWidth', 1.5,'Color',[0 0.4470 0.7410]); hold on
semilogy(time_hours, est_cov_lon, '--r', 'LineWidth', 1.5);
hold off
grid on
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Longitude Error [rad]', 'Interpreter', 'latex', 'FontSize', 12);
title('Longitude Error and 3-Sigma Bound (time evolution in MOONLANDER\_TOPO)', 'Interpreter', 'latex', 'FontSize', 10);
legend({'Longitude error', '$3\sigma_{\mathrm{lon}}$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');

fprintf('3.4)----------------------SC state and lunar lander coordinates---------------------\n')
fprintf('----------------------using position vector and range measurements------------------\n')

% Display lunar lander absolute state at initial time
fprintf('The estimated state at initial time is:\n');
disp(E_UKF_4(:, 1))
fprintf('\n')
fprintf('The estimated state at final time is:\n');
disp(E_UKF_4(:, end))
fprintf('With related covariance matrix:\n');
for i = 1:size(P_UKF_4, 1)
    fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P_UKF_4(i, :, end));
end
fprintf('\n');
fprintf('The more refined estimate of lander latitudinal coordinates are:\n')
fprintf('-latitude: LAT = %.4f [deg]\n', E_UKF_4(7, end)*cspice_dpr)
fprintf('-longitude: LON = %.4f [deg]\n', E_UKF_4(8, end)*cspice_dpr)
fprintf('-altitude: ALT = %f [m]\n', 0)
fprintf('\n\n\n')

%% Functions

%--------------------------------------------------------------------------
function dSdt = TBP(~, S, mu)
%
% DESCRIPTION:
% This function computes the time derivative of the state vector for the
% Two-Body Problem (TBP), where the motion of a body is governed solely by
% the gravitational attraction of a central body.
%
% PROTOTYPE
% dSdt = TBP(~, S, mu)
% Inputs:
%   ~  - [scalar] Time (not used since TBP is autonomous).
%   S  - [6x1] State vector [x; y; z; vx; vy; vz], where:
%          x, y, z   [km]     Position components in Cartesian coordinates.
%          vx, vy, vz [km/s]  Velocity components in Cartesian coordinates.
%   mu - [scalar] Gravitational parameter [km^3/s^2] of the central body.
%
% Output:
%   dSdt - [6x1] Time derivative of the state vector, where:
%          dx/dt, dy/dt, dz/dt are the velocity components [km/s],
%          dvx/dt, dvy/dt, dvz/dt are the acceleration components [km/s^2].

    % Compute the inverse cube of the distance from the central body
    rr = S(1:3);

    % Compute the square distance and distance:
    r2 = dot(rr, rr);
    r = sqrt(r2);

    % Assemble the time derivative of the state vector
    dSdt = [S(4:6); -mu*rr/(r2*r)];
end

%--------------------------------------------------------------------------
function [El_land, vw, xx_orb_MCI] = visibility_windows(S0_orb, tspan, par)
%
% DESCRIPTION:
% This function determines the time intervals during which an orbiter is
% visible from a lunar lander, based on their relative positions and a
% specified minimum elevation angle (mask angle). It propagates the
% orbiter's trajectory, calculates its position relative to the lander in
% a topocentric reference frame, and identifies visibility times.
%
% PROTOTYPE:
% [El_land, vw, xx_orb_MCI] = visibility_windows(S0_orb, tspan, par)
%
% Inputs:
%   - S0_orb: [6x1] Initial state vector of the orbiter in Moon-Centered Inertial (MCI) frame.
%             [x, y, z, vx, vy, vz], position in km, velocity in km/s.
%   - tspan: [1xN] Time vector for propagation (ET seconds).
%   - par: [struct] Structure containing parameters:
%           * mu_M: Gravitational parameter of the Moon [km^3/s^2].
%           * mask_angle: Minimum elevation angle for visibility [rad].
%           * LL0_land: [3x1] Latitudinal coordinates of the lander in IAU_MOON frame
%             (range [km], longitude [rad], latitude [rad]).
%
% Outputs:
%   - El_land: [1xM] Elevation coordinates of the orbiter in the lander’s 
%              topocentric reference frame during visibility.
%   - vw: [1xM] Vector of times (ET seconds) during which the orbiter is visible from the lander.
%   - xx_orb_MCI: [6xN] Orbiter state trajectory in the J2000 (MCI) frame.

    % Step 1: Extract parameters from the input structure
    mu = par.mu_M;               % Gravitational parameter of the Moon
    mask_angle = par.mask_angle; % Minimum elevation angle for visibility
    lat = par.LL0_land(2);       % Latitude of the station [rad]
    lon = par.LL0_land(3);       % Longitude of the station [rad]
    alt = 0;                     % Altitude of the station [km]

    % Compute the flattening:
    flat_MOON = (par.M_Radii(1) - par.M_Radii(3))/par.M_Radii(1);

    % Compute station position wrt Moon center (in MCMF or IAU_MOON)
    rv0_land_IAUMOON = cspice_pgrrec('MOON', lon, lat, alt, par.M_Radii(1), flat_MOON);

    % Step 2: Define number of time steps and numerical integration options
    mt = length(tspan); % Number of time steps
    options_ode = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3)]); 

    % Step 5: Propagate orbiter's trajectory in the MCI frame using TBP
    [~, xx_orb_MCI] = ode78(@(t, x) TBP(t, x, mu), tspan, S0_orb, options_ode);

    % Step 6: Compute the rotation matrix from IAU_MOON to J2000 at each time step
    ROT_IAUMOON2MCI = cspice_pxform('IAU_MOON', 'J2000', tspan);

    % Step 7: Rotate lander's position to J2000 frame
    rv0_land_MCI = pagemtimes(ROT_IAUMOON2MCI, reshape(repmat(rv0_land_IAUMOON, 1, mt), size(rv0_land_IAUMOON, 1), 1, mt));
    rv0_land_MCI = squeeze(rv0_land_MCI);

    % Step 8: Compute relative position of the orbiter w.r.t. the lander in MCI frame
    rv_orbland_MCI = xx_orb_MCI(:, 1:3)' - rv0_land_MCI;

    % Step 9: Compute rotation matrix from MCI to the lander's topocentric frame
    % 1) 9.1: Compute IAUMOON2TOPO and MCI2IAUMOON rotation matrix
    R_IAUMOON2MOONLANDER = cspice_eul2m(lat - pi, pi - lon, pi / 2, 2, 1, 2);
    R_MCI2IAUMOON = cspice_pxform('J2000', 'IAU_MOON', tspan);

    % 2) 9.2: Compute the rotation matrix from MCI to the lander's topocentric frame
    ROT_MCI2MOONLANDER = pagemtimes(repmat(R_IAUMOON2MOONLANDER, 1, 1, mt), R_MCI2IAUMOON);

    % Step 10: Rotate relative position to the lander's topocentric frame
    rv_orbland_topo = pagemtimes(ROT_MCI2MOONLANDER, reshape(rv_orbland_MCI, 3, 1, mt));
    rv_orbland_topo = squeeze(rv_orbland_topo);

    % Step 11: Convert relative position to latitudinal coordinates in topocentric frame
    [~, ~, Elevation] = cspice_reclat(rv_orbland_topo(1:3, :));

    % Step 12: Identify visibility times based on elevation angle
    visible_indices = Elevation >= mask_angle; % Logical array of visibility times

    % Step 13: Process visibility results
    if any(visible_indices)
        % Extract visibility times and coordinates
        vw = tspan(visible_indices);          % Times of visibility
        El_land = Elevation(visible_indices); % Visible elevations
    else
        % Handle case with no visibility windows
        fprintf('There is no visibility window. \n');
    end
end

%--------------------------------------------------------------------------
function [E, P] = UKF(S0, P0, tspan, meas, par, case_UKF)
%
% DESCRIPTION
% Apply the Unscented Kalman Filter to the measurements (meas), over the
% time span (tspan), starting from the initial mean state S0, and initial
% covariance matrix P0. 
%
% PROTOTYPE:
% [E, P] = UKF(S0, P0, tspan, meas, par, case_UKF)
%
% INPUTS:
%   - S0: Initial state vector of the orbiter [nx1].
%   - P0: Initial covariance matrix of the state [nxn].
%   - tspan: Time vector for propagation (ET seconds) [1xM].
%   - meas: Noisy measurements of orbiter's position in the lander's 
%           topocentric reference frame [3xM].
%   - par: Structure containing filter parameters:
%       - par.mu_M: Gravitational parameter of the Moon.
%       - par.alpha: Scaling parameter for sigma point generation.
%       - par.beta: Beta parameter for UKF weighting.
%       - par.noise_cov: Covariance matrix of measurement noise [3x3].
%   - case_UKF: Defines the update scheme to be used in the filter.
%
% OUTPUTS:
%   - E: Estimated state at each time step [nxM].
%   - P: State covariance matrix at each time step [nxnxM].

    % Extract necessary parameters.
    mu = par.mu_M;             % Gravitational parameter of the Moon.
    alpha = par.alpha;         % Scaling parameter for sigma points.
    beta = par.beta;           % Beta parameter for weight adjustments.

    % Initialization and memory allocation.
    n = length(S0);          % State vector dimension.
    mt = length(tspan);      % Number of time steps.
    E = zeros(n, mt);        % Estimated states.
    P = zeros(n, n, mt);     % Covariance matrices.
    E(:, 1) = S0;            % Initialize state.
    P(:, :, 1) = P0;         % Initialize covariance.

    % Unscented transform parameters.
    k = 0;                           % Secondary scaling parameter.
    lambda = alpha^2 * (n + k) - n;  % Scaling factor.
    n_lambda = n + lambda;           % Adjusted dimension.

    % Compute weights for sigma points.
    Wm = [lambda / n_lambda; repmat(1/(2*n_lambda), 2*n, 1)]; % Mean weights.
    Wc = Wm;                                                  % Covariance weights.
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);                     % Adjust the first weight.

    % ODE solver settings.
    options_ode = odeset('RelTol', 1e-12, 'AbsTol', 1e-13);

    % UKF loop for each time step.
    for i = 2:mt
        % Generate sigma points from the previous state and covariance.
        P_sqrt = sqrtm(n_lambda * P(:, :, i-1));
        Chi_old = [E(:, i-1), E(:, i-1) + P_sqrt, E(:, i-1) - P_sqrt];

        % Propagate each sigma point through the system dynamics.
        Chi_new = zeros(size(Chi_old));
        for j = 1:size(Chi_old, 2)
            [~, xx] = ode78(@(t, x) TBP(t, x, mu), [tspan(i-1), tspan(i)], Chi_old(1:6, j), options_ode);
            Chi_new(1:6, j) = xx(end, :)';

            if strcmp(case_UKF, '2')
                Chi_new(7:8, j) = Chi_old(7:8, j);
            end
        end

        % Calculate the predicted mean state:
        x_k_minus = Chi_new * Wm;             

        % Measurement update.
        switch case_UKF

            case '1'
                % Deviation of sigma points
                diffx = Chi_new - x_k_minus;  

                % Predicted covariance of sigma points:
                P_k_minus = diffx * diag(Wc) * diffx';    

                % Directly use the propagated sigma points for measurements.
                y_k = x_k_minus(1:3);              % Predicted measurement mean.
                diffy = Chi_new(1:3, :) - y_k;     % Measurement deviations.
                Pyy_k = diffy * diag(Wc) * diffy' + par.noise_cov(1:3, 1:3); % Measurement covariance.

            case '2'
                % Measurement update
                diffx = [Chi_new(1:6, :) - x_k_minus(1:6); angdiff(repmat(x_k_minus(7:8), 1, size(Chi_new, 2)), Chi_new(7:8, :))]; % Deviation of sigma points
        
                % Predicted covariance of sigma points:
                P_k_minus = diffx * diag(Wc) * diffx';    
        
                % Time evolution of the relative range between the orbiter and the lander:
                rv0_land_IAUMOON = cspice_latrec(repmat(par.M_Radii(1), 1, size(Chi_new, 2)), Chi_new(8, :), Chi_new(7, :)); 
                
                % Compute the rotation matrix from IAU_MOON to J2000 at each time step
                ROT_IAUMOON2MCI = cspice_pxform('IAU_MOON', 'J2000', tspan(i));
        
                % Rotate the vector position of lander to MCI:
                rv0_land_MCI = pagemtimes(ROT_IAUMOON2MCI, reshape(rv0_land_IAUMOON, 3, 1, size(rv0_land_IAUMOON, 2)));
                rv0_land_MCI = squeeze(rv0_land_MCI);
        
                % Calculate the relative position in MCI:
                rv0_orbland_MCI = Chi_new(1:3, :) - rv0_land_MCI;
        
                % Transform sigma points into measurement space.
                Ups_k = [Chi_new(1:3, :); vecnorm(rv0_orbland_MCI)];
                
                % Compute mean and deviations in measurement space.
                y_k = Ups_k * Wm;                                  % Predicted measurement mean.
                diffy = Ups_k - y_k;                               % Measurement deviations.
                Pyy_k = diffy * diag(Wc) * diffy' + par.noise_cov; % Measurement covariance.

            otherwise
                error('Invalid case for UKF.');
        end

        % Compute the cross-covariance and Kalman gain.
        Pxy_k = diffx * diag(Wc) * diffy';
        K_k = Pxy_k / Pyy_k;

        % Update state and covariance with the measurement.
        E(:, i) = x_k_minus + K_k * (meas(i, :)' - y_k);    % Posterior state
        P(:, :, i) = P_k_minus - K_k * Pyy_k * K_k';        % Posterior covariance
    end
end
