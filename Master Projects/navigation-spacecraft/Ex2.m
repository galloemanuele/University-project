% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2, Exercise 2
% Author: Emanuele Gallo

% Prepare the workspace:
clc; clearvars; cspice_kclear; close all

% Add path to sgp4:
addpath('sgp4');

% Load kernels:
cspice_furnsh('assignment02.tm');

% Set default seed:
rng('default')

% Some parameters useful for sgp4:
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Impose long format to better visualize results:
format long g

% Report the numerical values contained in Table 2, that is list of 
% stations and their corresponding measurement frequencies:
station.name = {'KOUROU', 'TROLL', 'SVALBARD'};
station.frequencies = [60, 30, 60];                                                 % [s]: Measurement frequencies in seconds for each station
station.mask_angle = [6*cspice_rpd, 0, 8*cspice_rpd];                               % [rad]: Mask angle for the stations
station.cov = diag([(125e-3*cspice_rpd).^2, (125e-3*cspice_rpd).^2, 0.01.^2]);      % [rad, rad, km]: Measurement noise
station.cpp = [30e3, 35e3, 35e3];                                                   % [€]: Cost per pass from each ground station
station.latitude = [5.25144, -72.011977, 78.229772];                                % [deg]: Stations latitude
station.longitude = [-52.80466, 2.536103, 15.407786];                               % [deg]: Station longitude


% Set option for propagations:
options_ode = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3)]);  % Options for ode78 solver


%% Exercise 2.1
% Define the spacecraft name and ID:
spacecraftName = 'SMOS';
spacecraftID = 36036;

% Corresponding strings in the TLE:
longstr1 = '1 36036U 09059A   24323.76060260  .00000600  00000-0  20543-3 0  9995';
longstr2 = '2 36036  98.4396 148.4689 0001262  95.1025 265.0307 14.39727995790658';

% Initialise the satrec structure:
satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

% Extract the gravitational constant from the satrec:
mu_E = satrec.mu;

% Convert Julian Date to calendar date and time:
[year,mon,day,hr,minutes,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);

% Convert date and time to SPICE ephemeris time (ET):
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,minutes,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% Print the satellite ID and the reference epoch:
fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s \n', sat_epoch_str);

% Update to the reference date the parameter for nutation correction:
arcsec2rad = pi/(180*3600);   % [rad]: Constant for arcseconds to radians conversions
ddpsi = -0.115178*arcsec2rad; % [rad]
ddeps = -0.007395*arcsec2rad; % [rad]

% Evaluate the TLE in True Equator Mean Equinox (TEME) reference frame:
[satrec, rteme_ref, vteme_ref] = sgp4(satrec, 0.0);

% Convert from TEME to ECI: 
ttt_ref = cspice_unitim(sat_epoch_et, 'ET', 'TDT')/cspice_jyear()/100;                              % Centuries from TDT 2000 January 1 00:00:00.000 
[reci_ref, veci_ref, aeci_ref] = teme2eci(rteme_ref, vteme_ref, [0; 0; 0], ttt_ref, ddpsi, ddeps);  % Perform the conversion
Seci_ref = [reci_ref; veci_ref];

% Compute the osculating elements:
osc_el_ref = cspice_oscelt([rteme_ref; vteme_ref], sat_epoch_et, mu_E);
a_ref = osc_el_ref(1)/(1-osc_el_ref(2));                            % [km]: Semi-major axis
e_ref = osc_el_ref(2);                                          % [-]: Eccentricity
i_ref = osc_el_ref(3);                                          % [rad]: Inclination
raan_ref = osc_el_ref(4);                                       % [rad]: Longitude of ascending node 
om_ref = osc_el_ref(5);                                         % [rad]: Argument of periapsis
theta_ref = osc_el_ref(6);                                      % [rad]: True anomaly 
fprintf('Keplerian Elements:\n');
fprintf('Semi-major axis (a)       = %.3f km\n', a_ref);
fprintf('Eccentricity (e)         = %.6f\n', e_ref);
fprintf('Inclination (i)          = %.6f rad (%.3f deg)\n', i_ref, i_ref*cspice_dpr());
fprintf('Longitude of ascending node (RAAN) = %.6f rad (%.3f deg)\n', raan_ref, raan_ref*cspice_dpr());
fprintf('Argument of periapsis (ω) = %.6f rad (%.3f deg)\n', om_ref, om_ref*cspice_dpr());
fprintf('True anomaly (θ)         = %.6f rad (%.3f deg)\n\n\n', theta_ref, theta_ref*cspice_dpr());

% Initial time:
t0_str = '2024-11-18T20:30:00.000';
et0 = cspice_str2et(t0_str);

% Final time:
tf_str = '2024-11-18T22:15:00.000';
etf = cspice_str2et(tf_str);

% Equally spaced vector, with appropriate step for each station:
npoints = [round((etf-et0)/station.frequencies(1))+1, round((etf-et0)/station.frequencies(2))+1, round((etf-et0)/station.frequencies(3))+1];
et_vect = {linspace(et0, etf, npoints(1)); linspace(et0, etf, npoints(2)); linspace(et0, etf, npoints(3))};

% Propagate the initial state (given in reference time) up to the initial time:
[~, xx] = ode78(@(t,x) TBP(t, x, mu_E), [sat_epoch_et et0], Seci_ref, options_ode);   % Propagate
Seci = xx(end, :)';

% Loop through each station to compute its visibility window, store results, and plot
l = length(station.name);   % Length of the time vector
fprintf('--------------------------Visibility windows--------------------\n');
for i = 1:l
    stationName = station.name{i};   % Station name
    et_vect_i = et_vect{i};          % Time vector for the current station
    
    % Compute visibility window for the stations and the relative azimuth,
    % elevation and range measurements:
    [Az, El, R, vw] = visibility_window(station, i, et_vect_i, Seci, mu_E, 'TBP');  
    results.(stationName).vw = vw;      % Allocate visibility window

    % Convert start and finish visibility ephemeris times to TDB:
    visibility_start = cspice_timout(vw(1), 'YYYY-MM-DD HR:MN:SC.###');
    visibility_end = cspice_timout(vw(end), 'YYYY-MM-DD HR:MN:SC.###');
    results.(stationName).visibility_start = visibility_start;
    results.(stationName).visibility_end = visibility_end;
    fprintf('%s Visibility Window: \n Start = %s (UTC time), \n End = %s (UTC time)\n', stationName, visibility_start, visibility_end);

    % Convert Azimuth and Elevation from radians to degrees and store all the values:
    results.(stationName).Az = Az;
    results.(stationName).El = El;
    results.(stationName).R = R;
    results.(stationName).Az_deg = Az*cspice_dpr();
    results.(stationName).El_deg = El*cspice_dpr();

    % Plot Azimuth vs Elevation for the current station
    figure; 
    plot(results.(stationName).Az_deg, results.(stationName).El_deg, ...
         'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'HandleVisibility', 'off');
    hold on;
    scatter(results.(stationName).Az_deg, results.(stationName).El_deg, 30, 'o', ...
            'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', ...
            'DisplayName', 'TBP position');
    yline(station.mask_angle(i) * cspice_dpr, 'LineWidth', 1, 'Color', 'r', ...
          'LineStyle', '--', 'Label', '$\mathrm{Minimum\ elevation}$', ...
          'LabelHorizontalAlignment', 'center', ...
          'LabelVerticalAlignment', 'middle', ...
          'Interpreter', 'latex');
    xlabel('$\mathrm{Azimuth [deg]}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Elevation [deg]}$', 'Interpreter', 'latex', 'FontSize', 14);
    title([stationName, ' Station: Az vs El'], 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    xlim([min(results.(stationName).Az_deg) - 5, max(results.(stationName).Az_deg) + 5]);
    ylim([min(results.(stationName).El_deg) - 2, max(results.(stationName).El_deg) + 5]);
end
hold off;

%% Exercise 2.2
% Initialize a structure to hold results for each station
results_m = struct();

for i = 1:l
    stationName = station.name{i};                  % Station name
    et_vect_visible = results.(stationName).vw;     % Current visibility window analyzed

    % Run the measurement function only in the visibility windows:
    [Az_m, El_m, R_m, vw_new] = measurements(station, i, et_vect_visible, satrec, ddpsi, ddeps, 'measurement');

    % Allocate the measurement visibility window:
    results_m.(stationName).vw = vw_new;      % Allocate visibility window

    % Check if the visibility window changes with the measurements:
    if ~isequal(size(vw_new), size(results.(stationName).vw))
        disp('The predicted visibility windows and the one from the measurement do not coincide');
        fprintf('For station %s \n', stationName)
    end

    % Convert Azimuth and Elevation from radians to degrees and store all the values:
    results_m.(stationName).Az = Az_m;
    results_m.(stationName).El = El_m;
    results_m.(stationName).R = R_m;
    results_m.(stationName).Az_deg = Az_m * cspice_dpr();
    results_m.(stationName).El_deg = El_m * cspice_dpr();

    % Plot Azimuth vs Elevation for the current station
    figure;
    plot(results.(stationName).Az_deg, results.(stationName).El_deg, ...
         'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'HandleVisibility', 'off');
    hold on;
    scatter(results.(stationName).Az_deg, results.(stationName).El_deg, 30, 'o', ...
            'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', 'none', ...
            'DisplayName', 'TBP position');
    plot(results_m.(stationName).Az_deg, results_m.(stationName).El_deg, ...
         'LineWidth', 1.5, 'Color', [0 1 1], 'LineStyle', '--', 'HandleVisibility', 'off');
    scatter(results_m.(stationName).Az_deg, results_m.(stationName).El_deg, 30, 'o', ...
            'MarkerEdgeColor', [0 1 1], 'MarkerFaceColor', 'none', ...
            'DisplayName', 'Measurement');
    xlabel('$\mathrm{Azimuth [deg]}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Elevation [deg]}$', 'Interpreter', 'latex', 'FontSize', 14);
    yline(station.mask_angle(i) * cspice_dpr, 'LineWidth', 1.5, 'Color', 'r', ...
          'LineStyle', '--', 'Label', '$\mathrm{Minimum\ elevation}$', ...
          'LabelHorizontalAlignment', 'center', ...
          'LabelVerticalAlignment', 'middle', ...
          'HandleVisibility', 'off', ...
          'Interpreter', 'latex');
    title([stationName, ' Station: Az vs El'], 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    legend('show', 'Interpreter', 'latex', 'Fontsize', 14)
    xlim([min(results_m.(stationName).Az_deg) - 5, max(results_m.(stationName).Az_deg) + 5]);
    ylim([min(results.(stationName).El_deg) - 2, max(results_m.(stationName).El_deg) + 5]);
end
hold off;



%% Exercise 2.3

% Define useful quantities to start the least square solution:
Wm = inv(sqrtm(station.cov));    % Weights matrix
x0_guess = Seci;                 % Initial guess
max_size = max([size(results_m.KOUROU.Az, 2), size(results_m.TROLL.Az, 2), size(results_m.SVALBARD.Az, 2)]);
diff_1 = max_size - size(results_m.KOUROU.Az, 2);
diff_2 = max_size - size(results_m.TROLL.Az, 2);
diff_3 = max_size - size(results_m.SVALBARD.Az, 2);

% Allocate the measurements:
meas_real(:,:,1) = [results_m.KOUROU.Az, zeros(1, diff_1); 
                    results_m.KOUROU.El, zeros(1, diff_1); 
                    results_m.KOUROU.R, zeros(1, diff_1)];
meas_real(:,:, 2) = [results_m.TROLL.Az, zeros(1, diff_2); 
                     results_m.TROLL.El, zeros(1, diff_2); 
                     results_m.TROLL.R, zeros(1, diff_2)];
meas_real(:,:,3) = [results_m.SVALBARD.Az, zeros(1, diff_3); 
                    results_m.SVALBARD.El, zeros(1, diff_3); 
                    results_m.SVALBARD.R, zeros(1, diff_3)];

% Allocate the time span used:
tspan_mat = [results_m.KOUROU.vw, zeros(1, diff_1); 
             results_m.TROLL.vw, zeros(1, diff_2); 
             results_m.SVALBARD.vw, zeros(1, diff_3)];


% sgp4 propagation from t_ref to t0:
tsince_23 = (et0 - sat_epoch_et) / 60.0;
[~, rteme_23, vteme_23] = sgp4(satrec, tsince_23);                                       % Propagate
ttt_et0 = cspice_unitim(et0, 'ET', 'TDT')/cspice_jyear()/100;                            % Centuries from TDT 2000 January 1 00:00:00.000 
[reci_23, veci_23, ~] = teme2eci(rteme_23, vteme_23, [0; 0; 0], ttt_et0, ddpsi, ddeps);  % Convert TEME to ECI
Seci_23 = [reci_23; veci_23];                                                            % State in ECI (at t0)
osc_el_0 = cspice_oscelt([rteme_23; vteme_23], et0, mu_E);                               % Convert the state in osculating elements
a_0 = osc_el_0(1)/(1-osc_el_0(2));                                                       % Exact semi-major axis at instant 0
i_0 = osc_el_0(3);                                                                       % Exact inclination at instant t0

% Options for the least square (minimu variance) navigation solution:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');     % Set lqnonlin options

% 2.3 case a)
idxa = 1;
fun_a = @(x) costfunction(x, tspan_mat, et0, Wm, station, idxa, meas_real(:,:,1), mu_E, 'TBP');                % Function to be minimised
[xsol_a, resnorm_a, residual_a, exitflag_a, ~, ~, jacobian_a] = lsqnonlin(fun_a, x0_guess, [], [], options);   % least square solution for case a
if exitflag_a>0 % post processing only if the convergence is reached
    [lsq.a.err_pos, lsq.a.err_vel, lsq.a.sqrt_trace_pos, lsq.a.sqrt_trace_vel, lsq.a.std_a, lsq.a.a, lsq.a.std_i, lsq.a.i, ~, ~] = ...
        post_processing(xsol_a, et0, jacobian_a, resnorm_a, residual_a, Seci_23, mu_E, ddpsi, ddeps, 'a');
else 
    fprintf('No solutions found for case a')
end


% 2.3 case b)
idxb = 1:3;
fun_b = @(x) costfunction(x, tspan_mat, et0, Wm, station, idxb, meas_real, mu_E, 'TBP');                      % Function to be minimised
[xsol_b, resnorm_b, residual_b, exitflag_b, ~, ~, jacobian_b] = lsqnonlin(fun_b, x0_guess, [], [], options);  % least square solution for case b
if exitflag_b>0 % post processing only if the convergence is reached
    [lsq.b.err_pos, lsq.b.err_vel, lsq.b.sqrt_trace_pos, lsq.b.sqrt_trace_vel, lsq.b.std_a, lsq.b.a, lsq.b.std_i, lsq.b.i, ~, ~] = ...
        post_processing(xsol_b, et0, jacobian_b, resnorm_b, residual_b, Seci_23, mu_E, ddpsi, ddeps, 'b');
else 
    fprintf('No solutions found for case b')
end

% 2.3 case c)
idxc = 1:3;
fun_c = @(x) costfunction(x, tspan_mat, et0, Wm, station, idxc, meas_real, mu_E, 'J2_TBP');                   % Function to be minimised
[xsol_c, resnorm_c, residual_c, exitflag_c, ~, ~, jacobian_c] = lsqnonlin(fun_c, x0_guess, [], [], options);  % least square solution for case c
if exitflag_c>0 % post processing only if the convergence is reached
    [lsq.c.err_pos, lsq.c.err_vel, lsq.c.sqrt_trace_pos, lsq.c.sqrt_trace_vel, lsq.c.std_a, lsq.c.a, lsq.c.std_i, lsq.c.i, ~, ~] = ...
        post_processing(xsol_c, et0, jacobian_c, resnorm_c, residual_c, Seci_23, mu_E, ddpsi, ddeps, 'c');
else 
    fprintf('No solutions found for case c')
end

% Long-time propagation:
tf_prop_str = '2024-11-19T08:30:00.000';                                  % 12 hours final propagation in str format
et_f_prop = cspice_str2et(tf_prop_str);                                   % Convert it into ephemeris time
et_vect_prop = linspace(et0, et_f_prop, 1000);                            % Create the vector
tsince_traj = (et_vect_prop - sat_epoch_et) / 60.0;                       % Vector of time for the sgp4 propagation
ttt_et_vect = cspice_unitim(tsince_traj, 'ET', 'TDT')/cspice_jyear()/100; % Centuries from TDT 2000 January 1 00:00:00.000 

% Properly allocate the vectors:
reci_traj = zeros(size(tsince_traj, 2), 3);   % Position vector
veci_traj = zeros(size(tsince_traj, 2), 3);   % Velocity vector
t_vect_date = cell(1, length(et_vect_prop));  % Allocate cell array of proper size

% Calculate the (assumed) exact trajectory:
for i = 1:size(tsince_traj, 2)
    [~, rteme_traj, vteme_traj] = sgp4(satrec, tsince_traj(i));                                                         % Propagate
    [reci_traj(i, :), veci_traj(i, :), ~] = teme2eci(rteme_traj, vteme_traj, [0; 0; 0], ttt_et_vect(i), ddpsi, ddeps);  % Convert TEME to ECI
   
end

% Propagate the different results according to the model they were obtained:
[~, xx_a] = ode78(@(t,x) TBP(t, x, mu_E), et_vect_prop, xsol_a, options_ode);      % Propagation for case a
[~, xx_b] = ode78(@(t,x) TBP(t, x, mu_E), et_vect_prop, xsol_b, options_ode);      % Propagation for case b
[~, xx_c] = ode78(@(t,x) J2_TBP(t, x, mu_E), et_vect_prop, xsol_c, options_ode);   % Propagation for case c

% Calculate the errors on the position for different cases:
err_pos_a = vecnorm(xx_a(:, 1:3) - reci_traj, 2, 2);
err_pos_b = vecnorm(xx_b(:, 1:3) - reci_traj, 2, 2);
err_pos_c = vecnorm(xx_c(:, 1:3) - reci_traj, 2, 2);

% Calculate the errors on the velocity for different cases:
err_vel_a = vecnorm(xx_a(:, 4:6) - veci_traj, 2, 2);
err_vel_b = vecnorm(xx_b(:, 4:6) - veci_traj, 2, 2);
err_vel_c = vecnorm(xx_c(:, 4:6) - veci_traj, 2, 2);

% Vector of time for plot:
t_plot = (et_vect_prop - et_vect_prop(1))/3600;

% Plot the errors on the position:
figure; % Create a new figure
plot(t_plot, err_pos_a, 'LineWidth', 2)
hold on
plot(t_plot, err_pos_b, 'LineWidth', 2)
plot(t_plot, err_pos_c, 'LineWidth', 2)
hold off
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Position Error [km]', 'Interpreter', 'latex', 'FontSize', 14)
title('Position residual time evolution for different models (ECI frame)', 'Interpreter', 'latex', 'FontSize', 16)
legend({'KOUROUS (TBP)', 'All (TBP)', 'All (PTBP)'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 14)
grid on
set(gca, 'FontSize', 12)
xlim([0 12])
hold off

% Plot the errors on the velocity:
figure;
plot(t_plot, err_vel_a, 'LineWidth', 2)
hold on
plot(t_plot, err_vel_b, 'LineWidth', 2)
plot(t_plot, err_vel_c, 'LineWidth', 2)
hold off
xlabel('Time from $t_0$ [h]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Velocity Error [km/s]', 'Interpreter', 'latex', 'FontSize', 14)
title('Velocity residual time evolution for different models (ECI frame)', 'Interpreter', 'latex', 'FontSize', 16)
legend({'KOUROUS (TBP)', 'All (TBP)', 'All (PTBP)'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 14)
grid on
set(gca, 'FontSize', 12)
xlim([0 12])



%% 2.4

% The possible combinations are:
combinations = [1, 2; 1, 3; 2, 3];

% Pre-allocate the variables:
std = zeros(size(combinations));                                           % Standard deviations matrix
P = zeros(size(x0_guess, 1), size(x0_guess, 1), size(combinations, 1));    % Covariance matrix
P_ai = zeros(2, 2, size(combinations, 1));                                 % Covariance matrix for inclination and semi-major axis
a4 = zeros(size(combinations, 1));                                         % Semi-major axis
i4 = zeros(size(combinations, 1));                                         % Inclination
xsol = zeros(size(x0_guess, 1), size(x0_guess, 2), size(combinations, 1)); % Mean of the state
gen_std = zeros(size(combinations, 1), 1);                                 % Combined standard deviations for each combination

% Cycle over all the combinations:
for i = 1:size(combinations, 1)
    % Select one possible combination, run the lsqnonlin, and individuate the best solution:
    idx = combinations(i, :);
    fun = @(x) costfunction(x, tspan_mat, et0, Wm, station, idx, meas_real, mu_E, 'J2_TBP');                                        % Function to be minimised
    [xsol(:,:,i), resnorm, residual, exitflag, ~, ~, jacobian] = lsqnonlin(fun, x0_guess, [], [], options);                         % least square solution for case b

    if exitflag>0 % Post-processing if the algorithm converges
        [~, ~, ~, ~, std(i, 1), a4(i), std(i, 2), i4(i), P(:,:,i), P_ai(:,:,i)] = post_processing(xsol(:,:,i), et0, jacobian, resnorm, residual, Seci_23, mu_E, ddpsi, ddeps, 'no_print');   % Post-processing results elaboration
    else 
        error('No solution found')
    end

    % Compute the generalized variance of the combination:
    gen_std(i) = det(P_ai(:,:,i));
end

% Select the best combination based on the minimum combined uncertainty
[min_combined_std, min_comb] = min(gen_std);

% Display the results:
fprintf('\n\n---------------------------------Trade-off analysis------------------------------------\n');
fprintf('The best solution within the budget and the mission constraint is obtained passing from: \n');
fprintf('Station 1: %s\n', station.name{combinations(min_comb, 1)});
fprintf('Station 2: %s\n', station.name{combinations(min_comb, 2)});
fprintf('The navigation solution is: \n');
disp(xsol(:,:, min_comb));
fprintf('With associated covariance matrix: \n');
for i = 1:size(P, 1)
    fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P(i, :, min_comb));
end
fprintf('Standard deviation of Semi-major Axis (a): %.4e [km]\n', std(min_comb, 1));
fprintf('Standard deviation of Inclination (i): %.4e [deg]\n', std(min_comb, 2)*cspice_dpr);


%% 2.5 

% Set the one week propagation time:
tf5_str = '2024-11-21T20:30:00.000';
etf5 = cspice_str2et(tf5_str);

% Equally spaced vector, with appropriate step for each station:
npoints5 = [round((etf5-et0)/station.frequencies(1))+1, round((etf5-et0)/station.frequencies(2))+1, round((etf5-et0)/station.frequencies(3))+1];
et_vect5 = {linspace(et0, etf5, npoints5(1)); linspace(et0, etf5, npoints5(2)); linspace(et0, etf5, npoints5(3))};

% Loop through each station to compute its visibility window, store results, and plot
l = length(station.name);   % Length of the time vector
figure('Name', 'Ground Station Elevation vs Time', 'NumberTitle', 'off'); % Create a named figure

for i = 1:l
    stationName = station.name{i};   % Station name
    et_vect_i5 = et_vect5{i};        % Time vector for the current station

    % Compute visibility window for the stations and the relative azimuth,
    % elevation, and range measurements:
    [Az5, El5, R5, vw5] = visibility_window(station, i, et_vect_i5, Seci_23, mu_E, 'J2_TBP');
    
    % Allocate the visibility windows for this case:
    results5.(stationName).vw = vw5;

    % Convert Azimuth and Elevation from radians to degrees and store all the values:
    results5.(stationName).Az = Az5;
    results5.(stationName).El = El5;
    results5.(stationName).R = R5;
    results5.(stationName).Az_deg = Az5*cspice_dpr();
    results5.(stationName).El_deg = El5*cspice_dpr();

    % % Convert null values to nan for plotting enhancement reasons:
    El5(El5==0) = nan;

    % Plot Elevation vs Time for the current station
    subplot(3, 1, i);
    plot((et_vect_i5 - et_vect_i5(1))/(3600*24), El5*cspice_dpr, ...
         'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
         'DisplayName', sprintf('Elevation for station %s', stationName));
    hold on;
    yline(station.mask_angle(i) * cspice_dpr, 'LineWidth', 1, 'Color', 'r', ...
          'LineStyle', '--', 'LineWidth', 1.5);
    ylabel('El [deg]', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Time elapsed from $t_0$ [days]', 'Interpreter', 'latex', 'FontSize', 14);
    title(sprintf('%s Station', stationName), 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    ylim([-4, max(results5.(stationName).El_deg) + 5]);
    xlim([0, max(et_vect_i5 - et_vect_i5(1))/cspice_spd]);
end

% Add the overall title for all subplots
% sgtitle('Elevation vs Time for Ground Stations', 'Interpreter', 'latex', 'FontSize', 20);
hold off;
% From the previous plot, it can be clearly seen that the repeatibility of
% the measurement is achievable selecting at least two ground stations.
% From this plot, it is also noticeable that the KOUROU ground station
% cannot be selected as prime ground station, but maybe only as backup one.


%% Functions

%-----------------------------------------------------------------------
function dSdt = TBP(~, S, mu)
% TBP: Two-Body Problem (TBP) dynamics model.
%
% Description:
% This function computes the time derivative of the state vector for the
% Two-Body Problem (TBP), where the motion of a body is governed solely by
% the gravitational attraction of a central body.
%
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

%-----------------------------------------------------------------------
function dSdt = J2_TBP(et, S, mu)
% J2_TBP - J2-perturbed Two-Body Problem.
%
% Inputs:
%   et  - Ephemeris time.
%   S   - State vector [x; y; z; vx; vy; vz] (position in km, velocity in km/s).
%   mu  - Gravitational parameter (km^3/s^2) of the central body.
%
% Output:
%   dSdt - Time derivative of the state vector [dx/dt; dy/dt; dz/dt; dvx/dt; dvy/dt; dvz/dt].
    
    % Allocate the space:
    dSdt = zeros(6, 1);

    % Compute the inverse cube of the distance from the central body
    rr = S(1:3);

    % Compute the square distance and distance:
    r2 = dot(rr, rr);
    r = sqrt(r2);

    % Assemble the time derivative of the state vector
    dSdt(1:3) = S(4:6); 
    a_grav = -mu*rr/(r2*r);

    % Rotation matrix from ECI to ECEF:
    ROT_ECI2ECEF = cspice_pxform('J2000', 'ITRF93', et);

    % Position in ECEF:
    rr_ECEF = ROT_ECI2ECEF * rr;

    % Compute the square distance and distance:
    r2_ECEF = dot(rr, rr);
    r_ECEF = sqrt(r2_ECEF);

    % Impose the J2 value:
    J2 = 0.0010826269;

    % Compute the radius of the Earth:
    E_Radii = cspice_bodvrd('EARTH', 'RADII', 3);
    Re = E_Radii(1);

    % Compute the acceleration in ECEF:
    a_J2_ECEF = 3/2*mu*J2*rr_ECEF/(r2_ECEF*r_ECEF)*(Re/r_ECEF).^2.*(5*(rr_ECEF(3)/r_ECEF).^2-[1;1;3]);

    % Convert it to ECI:
    ROT_ECEF2ECI = cspice_pxform('ITRF93','J2000', et); % Rotation matrix from ECI to ECEF
    a_J2_ECI = ROT_ECEF2ECI * a_J2_ECEF;                % Convert acceleration to ECI

    % Add J2 acceleration to keplerian one:
    dSdt(4:6) = a_grav + a_J2_ECI;
end

%-----------------------------------------------------------------------
function [Az, El, R, vw] = visibility_window(station, idx, et_vect, Seci, mu, propagation_type)
    % VISIBILITY_WINDOW - Computes visibility parameters (azimuth, 
    % elevation, and range) of a spacecraft relative to a ground station 
    % for a given time vector using a vectorized approach.
    %
    % Inputs:
    %   stationName - String containing the name of the ground station
    %   et_vect     - Vector of ephemeris times (ET) for which visibility is computed.
    %   Seci        - State vector of the spacecraft in the ECI frame (position and velocity).
    %
    % Outputs:
    %   Az          - Azimuth of the spacecraft relative to the station (radians).
    %   El          - Elevation of the spacecraft relative to the station (radians).
    %   R           - Range from the station to the spacecraft (kilometers).

    % Station-specific details:
    stationName = station.name{idx};      % Extract the station name
    mask_angle = station.mask_angle(idx); % Extract the station mask angle
    topoFrame = [stationName, '_TOPO'];   % Topocentric frame name based on stationName

    % Number of time points:
    mt = length(et_vect);

    % Step 0: Propagate the spacecraft state:
    options_ode = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3)]);  % Options for ode78 solver
    switch propagation_type
        case 'TBP'
            [~, xx] = ode78(@(t,x) TBP(t, x, mu), et_vect, Seci, options_ode); 
        case 'J2_TBP'
            [~, xx] = ode78(@(t,x) J2_TBP(t, x, mu), et_vect, Seci, options_ode);
        otherwise 
            error('The propagation type chosen is not correct')
    end

    % Step 1: Compute station positions and velocities in ECI frame
    rv_station_eci = cspice_spkezr(stationName, et_vect, 'J2000', 'NONE', 'EARTH');

    % Step 2: Compute station-to-satellite vectors in the ECI frame
    rv_station_sat_eci = xx' - rv_station_eci; % 6xmt matrix

    % Step 3: Perform coordinate transformations from ECI to topocentric frames
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et_vect); 

    % Step 4: Transform states from ECI to topocentric frame
    rv_station_sat_topo = pagemtimes(ROT_ECI2TOPO, reshape(rv_station_sat_eci, 6, 1, mt)); 
    rv_station_sat_topo = squeeze(rv_station_sat_topo);                                    

    % Step 5: Compute range, elevation, and azimuth
    [Range, Azimuth, Elevation] = cspice_reclat(rv_station_sat_topo(1:3, :));

    % Step 6: Apply visibility mask angle:
    visible_indices = Elevation >= mask_angle; % Logical index for visibility

    if any(visible_indices)
        switch propagation_type
            case 'TBP'
                % Extract visibility window in ET seconds
                vw = et_vect(visible_indices);
        
                % Filter azimuth, elevation, and range for visible time points
                Az = Azimuth(visible_indices);
                El = Elevation(visible_indices);
                R = Range(visible_indices);

            case 'J2_TBP'
                % Initialize output arrays with NaN
                Az = zeros(size(et_vect));
                El = zeros(size(et_vect));
                R = zeros(size(et_vect));
                vw = zeros(size(et_vect));

                % Assign values only to visible_indices
                Az(visible_indices) = Azimuth(visible_indices);
                El(visible_indices) = Elevation(visible_indices);
                R(visible_indices) = Range(visible_indices);
                vw(visible_indices) = et_vect(visible_indices);
            otherwise 
                error('The propagation type chosen is not correct')
        end 
    else
        fprintf('No visibility detected for the %s station.\n', stationName);
        vw = [];
        Az = [];
        El = [];
        R = [];
    end
end

%-----------------------------------------------------------------------
function [Az_m, El_m, R_m, vw_m] = measurements(station, idx, et_vect_visible, satrec, ddpsi, ddeps, type)
%
% DESCRIPTION:
% Function to compute azimuth, elevation, and range measurements from a ground station to a satellite.
% Measurements are affected by noise based on the station's covariance matrix.
%
% PROTOTYPE:
% [Az_m, El_m, R_m, vw_m] = measurements(station, idx, et_vect_visible, satrec, ddpsi, ddeps, type)
%
% INPUTS:
%   station         - Struct containing station data (e.g., name, mask angle, covariance matrix).
%   idx             - Index of the station in the station struct.
%   et_vect_visible - Vector of ephemeris times when the satellite is visible.
%   satrec          - SGP4 satellite record structure.
%   ddpsi, ddeps    - Nutation corrections for Earth's orientation model.
%   type            - Choose between 'measurements' to better perform in
%                     measuremets field and 'long_term' to better perform 
%                     in long-term analysis.
%
% OUTPUTS:
%   Az_m            - Vector of azimuth measurements [rad].
%   El_m            - Vector of elevation measurements [rad].
%   R_m             - Vector of range measurements [km].

    % Extract station properties:
    stationName = station.name{idx};             % Name of the station
    meas_noise_cov = station.cov;                % Measurement noise covariance matrix
    mask_angle = station.mask_angle(idx);        % Minimum elevation angle for visibility [rad]
    topoFrame = [stationName, '_TOPO'];          % Topocentric reference frame name

    % Length of time vector:
    mt = length(et_vect_visible);

    % Compute satellite epoch in ET
    [year, mon, day, hr, minutes, sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
    sat_epoch_et = cspice_str2et(sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', ...
                        [year, mon, day, hr, minutes, sec])); 

    % Time since epoch in minutes
    tsince_vect = (et_vect_visible - sat_epoch_et) / 60.0;

    % Step 1: Propagate satellite states using SGP4 and convert TEME to ECI
    rteme = NaN(3, length(tsince_vect));
    vteme = NaN(3, length(tsince_vect));
    ttt = cspice_unitim(et_vect_visible, 'ET', 'TDT') / cspice_jyear() / 100; % Julian centuries from J2000
    reci = NaN(3, length(tsince_vect));
    veci = NaN(3, length(tsince_vect));
    for i = 1:mt
        [~, rteme(:,i), vteme(:,i)] = sgp4(satrec, tsince_vect(i));                                     % Propagate
        [reci(:,i), veci(:,i), ~] = teme2eci(rteme(:,i), vteme(:,i), [0; 0; 0], ttt(i), ddpsi, ddeps);  % Convert TEME to ECI

    end

    % Step 2: Compute station states in ECI (vectorized SPICE call)
    rv_station_eci = cspice_spkezr(stationName, et_vect_visible, 'J2000', 'NONE', 'EARTH');

    % Step 3: Relative states in ECI
    rv_station_sat_eci = [reci; veci] - rv_station_eci;

    % Step 4: Transform states to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et_vect_visible); % 6x6xmt matrix
    rv_station_sat_topo = pagemtimes(ROT_ECI2TOPO, reshape(rv_station_sat_eci, 6, 1, mt));
    rv_station_sat_topo = squeeze(rv_station_sat_topo);

    % Step 5: Calculate range, azimuth and elevation
    [Range, Azimuth, Elevation] = cspice_reclat(rv_station_sat_topo(1:3, :));

    % Step 6: Add measurement noise
    rng default                                         % Impose always the same 
    meas_mean = [Azimuth; Elevation; Range]';           % Mean measurement matrix
    noisy_meas = mvnrnd(meas_mean, meas_noise_cov)';    % Noisy measurements (with meas_noise_cov covariance matrix)

    % Step 7: Apply visibility condition:
    valid_idx = noisy_meas(2, :) >= mask_angle;    % Spacecraft visibile if the eleveation is bigger than the visibility

    if any(valid_idx) % Only allocate variables if the visibility condition is respected

        switch type
            case 'measurement'
                R_m = noisy_meas(3, valid_idx);                % Only allocate range noisy measurements in the visibility field
                El_m = noisy_meas(2, valid_idx);               % Only allocate elevation noisy measurements in the visibility field
                Az_m = noisy_meas(1, valid_idx);               % Only allocate azimuth noisy measurements in the visibility field
                vw_m = et_vect_visible(valid_idx);             % Visibility window of the measurement case

             case 'long_term'
                 % Initialize output arrays with NaN
                 Az_m = zeros(size(et_vect_visible));
                 El_m = zeros(size(et_vect_visible));
                 R_m = zeros(size(et_vect_visible));
                 vw_m = zeros(size(et_vect_visible));

                 % Assign values only to visible_indices
                 Az_m(valid_idx) = Azimuth(valid_idx);
                 El_m(valid_idx) = Elevation(valid_idx);
                 R_m(valid_idx) = Range(valid_idx);
                 vw_m(valid_idx) = et_vect_visible(valid_idx);

            otherwise 
                error('Incorrect measurement type chosen')
        end

        else
        fprintf('No visibility detected for the %s station.\n', stationName);

        % Do not allocate any variable:
        R_m = [];
        El_m = [];
        Az_m = [];
        vw_m = [];
    end
end

%-----------------------------------------------------------------------
function res = costfunction(x0, tspan_mat, t0, Wm, station, idx, meas_real_mat, mu, propagation_type)
% 
% DESCRIPTION
% Function to compute the residual between predicted and real measurements in a 
% trajectory propagation problem under the influence of the Two-Body Problem (TBP).
% 
% PROTOTYPE:
% res = costfunction(x0, tspan_mat, t0, Wm, station, idx, meas_real_mat, mu, propagation_type)
%
% INPUTS:
%   x0            - Initial state vector [position; velocity].
%   tspan_mat     - Double array of time spans for propagation.
%   t0            - Initial time for propagation.
%   Wm            - Weight matrix for residual computation.
%   station       - Struct containing station information.
%   idx           - Number of stations (or iterations).
%   meas_real_mat - Double array of real measurements for each station.
%   mu            - Gravitational parameter [km^3/s^2].
%   propagation_type  - Choose the type of propagation between 'TBP',
%                       and 'J2_TBP'
% OUTPUTS:
%   res         - Residuals between predicted and real measurements, weighted by Wm.

    % Properly pre-allocate the residual variable:
    res = nan(size(meas_real_mat, 1), size(meas_real_mat(idx, :), 2), max(idx));

    for i = idx     
        % Allocate the time variables:
        tspan = tspan_mat(i, :); % Define the time vector
        tspan = tspan(tspan~=0); % Allocate only non null elements
        mt = length(tspan);      % Time vector length      

        % Allocate the measurement variables:
        meas_real = meas_real_mat(:,:,i);                    % Choose the right measurement
        rows_non_null = any(meas_real ~= 0, 2);              % Logical vector for non-null rows
        cols_non_null = any(meas_real ~= 0, 1);              % Logical vector for non-null columns
        meas_real = meas_real(rows_non_null, cols_non_null); % measurements with null elements filtered-out

        % Station-specific details:
        stationName = station.name{i};
        topoFrame = [stationName, '_TOPO']; % Topocentric frame name based on stationName.
        
        % Compute the predicted measurements using numerical integration
        options_ode = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3)]);      % ODE solver options (set to extemely precise, to have good measurements)
        switch propagation_type
            case 'TBP'
            [~, xx0] = ode78(@(t,x) TBP(t, x, mu), [t0, tspan(1)], x0, options_ode);              % Propagate up to visibility begin 
            [~, xx_pre] = ode78(@(t,x) TBP(t, x, mu), tspan, xx0(end, :), options_ode);           % Propagate state using TBP

            case 'J2_TBP'
            [~, xx0] = ode78(@(t,x) J2_TBP(t, x, mu), [t0, tspan(1)], x0, options_ode);             % Propagate up to visibility begin 
            [~, xx_pre] = ode78(@(t,x) J2_TBP(t, x, mu), tspan, xx0(end, :), options_ode);          % Propagate state using PTBP

        otherwise % Non valid case
            fprintf('You have selected an invalid propagation type')
        end
        
        % Step 1: Compute station positions and velocities in ECI frame
        rv_station_eci = cspice_spkezr(stationName, tspan, 'J2000', 'NONE', 'EARTH'); % 6xmt matrix
    
        % Step 2: Compute station-to-satellite vectors in the ECI frame
        rv_station_sat_eci = xx_pre' - rv_station_eci; % 6xmt matrix
    
        % Step 3: Perform coordinate transformations from ECI to topocentric frames
        ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, tspan); % 6x6xmt matrix
    
        % Step 4: Transform states from ECI to topocentric frame
        rv_station_sat_topo = pagemtimes(ROT_ECI2TOPO, reshape(rv_station_sat_eci, 6, 1, mt)); % 6x1xmt matrix
        rv_station_sat_topo = squeeze(rv_station_sat_topo);                                    % Convert back to 6xmt matrix
    
        % Step 5: Compute range, elevation, and azimuth
        [R, Az, El] = cspice_reclat(rv_station_sat_topo(1:3, :));
    
        % Step 6: Allocate the measurements:
        meas_pred = [Az; El; R];
    
        % Step 7: Compute the residuals:
        res(:, 1:mt, i) = Wm * [angdiff(meas_pred(1:2, :), meas_real(1:2, :)); meas_pred(3, :) - meas_real(3, :)];
    end

    % Post-processing:
    res = reshape(res, size(res, 1), size(res, 3)*size(res, 2)); % Reshape the residual

    % Allocate only the non-nan values:
    res1 = res(1, :); res1 = res1(~isnan(res1));
    res2 = res(2, :); res2 = res2(~isnan(res2));
    res3 = res(3, :); res3 = res3(~isnan(res3));
    res = [res1; res2; res3];
end

%-----------------------------------------------------------------------
function [err_pos, err_vel, sqrt_trace_pos, sqrt_trace_vel, std_a, a, std_i, i, P, P_ai] = post_processing(xsol, et, jacobian, resnorm, residual, Seci_23, mu, ddpsi, ddeps, case_analysed)
%
% DESCRIPTION
% This function post-processes the information exitign from lsqnonlin to 
% compute and displays the navigation solution metrics, making them
% available for reading and/or comparison. 
%
% PROTOTYPE:
% [err_pos, err_vel, sqrt_trace_pos, sqrt_trace_vel, std_a, a, std_i, i, P, P_ai] = ... 
% post_processing(xsol, et, jacobian, resnorm, residual, Seci_23, mu, ddpsi, ddeps, case_analysed)
%
% INPUTS:
%   xsol          - State solution vector [km, km/s]
%   jacobian      - Jacobian matrix (sparse)
%   resnorm       - Residual norm
%   residual      - Residual vector
%   Seci_23       - State vector dimensions for the residual calculation
%   mu_E          - Gravitational parameter of Earth [km^3/s^2]
%   case_analysed - Case identifier ('a', 'b', or 'c')
% 
% OUTPUTS:

    
    % Allocate properly the variables:
    reci_23 = Seci_23(1:3);
    veci_23 = Seci_23(4:6);

    % Calculate the errors:
    err_pos = norm(reci_23 - xsol(1:3));
    err_vel = norm(veci_23 - xsol(4:6));
    
    % Covariance computation:
    n = length(residual);
    m = length(Seci_23);
    Jac = full(jacobian);
    JacJ = Jac.' * Jac;
    P = resnorm / (n - m) * (JacJ \ eye(size(JacJ)));
    sqrt_trace_pos = sqrt(trace(P(1:3, 1:3)));
    sqrt_trace_vel = sqrt(trace(P(4:6, 4:6)));
    
    % State elements:
    x = xsol(1); y = xsol(2); z = xsol(3);
    vx = xsol(4); vy = xsol(5); vz = xsol(6);
    
    % Jacobian for Keplerian elements (a, i), it is computed symbolically:
    Jac_ai_S = [ (2*mu^2*x)/((x^2 + y^2 + z^2)^(3/2)*(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2),...
                (2*mu^2*y)/((x^2 + y^2 + z^2)^(3/2)*(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2),...
                (2*mu^2*z)/((x^2 + y^2 + z^2)^(3/2)*(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2), ...
                (2*mu*vx)/(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2, ...
                (2*mu*vy)/(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2, ...
                (2*mu*vz)/(vx^2 - (2*mu)/(x^2 + y^2 + z^2)^(1/2) + vy^2 + vz^2)^2;...
                ((vz*y - vy*z)*(z*vx^2 - vz*x*vx + z*vy^2 - vz*y*vy))/((1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 ...
                + (vz*x - vx*z)^2 + (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2)),...
                -((vz*x - vx*z)*(z*vx^2 - vz*x*vx + z*vy^2 - vz*y*vy))/((1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 + ...
                (vz*x - vx*z)^2 + (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2)),...
                -((2*vx*(vz*x - vx*z) + 2*vy*(vz*y - vy*z))*(vy*x - vx*y))/(2*(1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 +...
                (vz*x - vx*z)^2 + (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2)),...
                ((vz*y - vy*z)*(vz*x^2 - vx*z*x + vz*y^2 - vy*z*y))/((1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 +...
                (vz*x - vx*z)^2 + (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2)),...
                -((vz*x - vx*z)*(vz*x^2 - vx*z*x + vz*y^2 - vy*z*y))/((1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 +...
                (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2)), ((2*x*(vz*x - vx*z) +...
                2*y*(vz*y - vy*z))*(vy*x - vx*y))/(2*(1 - (vy*x - vx*y)^2/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 +...
                (vz*y - vy*z)^2))^(1/2)*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2))];

    % Convert the navigation solution mean state from eci to teme:
    ttt = cspice_unitim(et, 'ET', 'TDT')/cspice_jyear()/100;                            % Centuries from TDT 2000 January 1 00:00:00.000 
    [rteme, vteme, ~] = eci2teme(xsol(1:3), xsol(4:6), [0; 0; 0], ttt, ddpsi, ddeps);   % Convert from eci to teme

    % Properly allocate the sami-major axis and inclination:
    osc_el = cspice_oscelt([rteme; vteme], et, mu);                 % Convert navigation solution into osculating elements
    a  = osc_el(1)/(1-osc_el(2));                         % [km]: Semi-major axis
    i = osc_el(3);                                        % [rad]: Inclination

    % Compute the covariance matrix associated to the reduced keplerian
    % mean state:
    P_ai = Jac_ai_S * P * Jac_ai_S';
    
    % Extract standard deviations for a and i:
    std_a = sqrt(P_ai(1, 1));  % Standard deviation of semimajor axis
    std_i = sqrt(P_ai(2, 2));  % Standard deviation of inclination
    
    % Display results based on case:
    switch case_analysed
        case 'a'
            fprintf('\n----------------- Navigation Solution for KOUROU (TBP)----------------------\n');
            fprintf('Analysis performed using measurements from the KOUROU ground station.\n The results are referred to initial time t0 and ECI frame. \n');
            fprintf('The navigation solution mean is: \n');
            disp(xsol)
            fprintf('With associated covariance matrix: \n');
            for j = 1:size(P, 1)
                fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P(j, :));
            end
            fprintf('Position Error (norm): %.10f [km]\n', err_pos);
            fprintf('Velocity Error (norm): %.10f [km/s]\n', err_vel);
            fprintf('Square root of the trace of position covariance: %.10f [km]\n', sqrt_trace_pos);
            fprintf('Square root of the trace of velocity covariance: %.10f [km/s]\n', sqrt_trace_vel);
            fprintf('Standard deviation of Semimajor Axis (a): %.4e [km]\n', std_a);
            fprintf('Standard deviation of Inclination (i): %.4e [deg]\n', std_i*cspice_dpr);
        case 'b'
            fprintf('\n------------ Navigation Solution for all stations (TBP)---------------------\n');
            fprintf('Analysis performed using simulated measurements from three ground stations.\n The results are referred to initial time t0 and ECI frame. \n');
            fprintf('The navigation solution mean is \n');
            disp(xsol)
            fprintf('With associated covariance matrix: \n');
            for j = 1:size(P, 1)
                fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P(j, :));
            end
            fprintf('Position Error (norm): %.10f [km]\n', err_pos);
            fprintf('Velocity Error (norm): %.10f [km/s]\n', err_vel);
            fprintf('Square root of the trace of position covariance: %.10f [km]\n', sqrt_trace_pos);
            fprintf('Square root of the trace of velocity covariance: %.10f [km/s]\n', sqrt_trace_vel);
            fprintf('Standard deviation of Semimajor Axis (a): %.4e [km]\n', std_a);
            fprintf('Standard deviation of Inclination (i): %.4e [deg]\n', std_i * cspice_dpr);
        case 'c'
            fprintf('\n------------- Navigation Solution for all stations (PTBP)--------------------\n');
            fprintf('Analysis performed using simulated measurements from the three gound stations \n considering J2-perturbed spacecraft dynamics.\n The results are referred to initial time t0 and ECI frame. \n');
            fprintf('The navigation solution mean is: \n');
            disp(xsol)
            fprintf('With associated covariance matrix: \n');
            for j = 1:size(P, 1)
                fprintf('[%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e]\n', P(j, :));
            end
            fprintf('Position Error (norm): %.10f [km]\n', err_pos);
            fprintf('Velocity Error (norm): %.10f [km/s]\n', err_vel);
            fprintf('Square root of the trace of position covariance: %.10f [km]\n', sqrt_trace_pos);
            fprintf('Square root of the trace of velocity covariance: %.10f [km/s]\n', sqrt_trace_vel);
            fprintf('Standard deviation of Semimajor Axis (a): %.4e [km]\n', std_a);
            fprintf('Standard deviation of Inclination (i): %.4e [deg]\n', std_i * cspice_dpr);
        case 'no_print'

        otherwise 
            fprintf('Invalid case selected as variable "case_analysed"');
    end

end
