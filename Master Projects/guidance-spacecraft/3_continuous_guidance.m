% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #1, Exercise 3
% Author: Emanuele Gallo

% Prepare the workspace:
clc; clearvars; cspice_kclear; close all

% Impose long format to better visualize results:
format long g

% Set dimensional parameters:
par.hi = 800;                               % Altitude of departure orbit [km]
par.i_i = 0;                                % Initial inclination [rad]
par.hf = 1000;                              % Altitude of arrival orbit [km]
par.Delta_i = deg2rad(0.75);                % Inclination change [rad]
par.i_f = deg2rad(par.i_i + par.Delta_i);   % Final inclination [rad]
par.Re = 6378.1366;                         % Earth radius [km]
par.rho_i = par.hi + par.Re;                % Initial distance from Earth [km]
par.rho_f = par.hf + par.Re;                % Final distance from Earth [km]
par.mu = 398600.435;                        % Earth gravitational parameter [km3/s2]
par.rho_0 = 750 + par.Re;                   % Reference radius for debris flux [km] 
par.k1 = 1e-5;                              % Debris spatial density constant 1 [DU^-1]
par.k2 = 1e-4;                              % Debris spatial density constant 2 [DU^2]
par.m0 = 1000;                              % Initial mass [kg]
par.T_max = 3;                              % Maximum thrust [N]
par.Isp = 3120;                             % Specific impulse [s]
par.g0 = 9.81;                              % Sea level gravity acceleration [m/s^2]
par.DU = 7178.1366;                         % Distance Unit [km]
par.MU = par.m0;                            % Mass Unit [kg]
par.TU = sqrt(par.DU^3/par.mu);             % Time Unit [s]
par.VU = par.DU/par.TU;                     % Velocity Unit [km/s]
par.FU = ((par.MU*par.DU*1000)/par.TU^2);   % Force Unit [N]
par.AU = (par.DU*1000)/(par.TU)^2;          % Acceleration Unit

%% Ex 3.1 

% Spatial density function:
q = @(rho) par.k1./(par.k2 + ((rho-par.rho_0)/par.DU).^2);

% Calculate the vector of distances from Earth:
rho_min = par.hi - 100 + par.Re;                       % Minimum distance from Earth [km]
rho_max = par.hf + 100 + par.Re;                       % Maximum distance from Earth [km]
rho_vect = linspace(rho_min, rho_max, 1000);           % Vector

% Calculate the distribution of debris spatial density:
q_rho = q(rho_vect);

% Find the maximum value of the density and the related value of distance
% from Earth:
[q_max,i] = max(q_rho);     % Maximum value of spatial debris density
rho_q_max = rho_vect(i);    % Dinstance from Earth for the maimum value of the spatial debris function

% Plot the distribution:
figure
plot(rho_vect, q_rho, 'LineWidth', 1.5, 'Color', [0.1, 0.5, 0.8]);
hold on

% Highlight the maximum point:
plot(rho_q_max, q_max, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.9, 0.3, 0.3], 'MarkerEdgeColor', 'k'); % Marker at max
xline(rho_q_max, '--', 'Color', [0.8, 0.1, 0.1], 'LineWidth', 1.5, 'HandleVisibility', 'off');            % Vertical line at max

% Highlight the departure orbit:
plot(par.rho_i, q(par.rho_i), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.3, 0.9, 0.3], 'MarkerEdgeColor', 'k'); % Departure point
xline(par.rho_i, '--', 'Color', [0.1, 0.8, 0.1], 'LineWidth', 1.5, 'HandleVisibility','off');                    % Vertical line for departure

% Highlight the arrival orbit:
plot(par.rho_f, q(par.rho_f), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.3, 0.3, 0.9], 'MarkerEdgeColor', 'k'); % Arrival point
xline(par.rho_f, '--', 'Color', [0.1, 0.1, 0.8], 'LineWidth', 1.5, 'HandleVisibility','off');                    % Vertical line for arrival
title('Spatial debris densisyt as a funciton of distance from Earth centre');
ylim([0, q_max + q_max/10]);
ax = gca;
ax.XTick = unique([ax.XTick, par.rho_i, par.rho_f, rho_q_max]); % Ensure key points are on x-axis ticks
ax.XAxis.TickLength = [0.015, 0.015];   % Set tick length
ax.YAxis.TickLength = [0.015, 0.015];
xlabel('Distance from Earth centre [km]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Spatial debris density [$DU^{-3}$]', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
ax.FontSize = 12;                  % Increase font size for better readability
ax.GridLineStyle = '--';           % Use dashed grid lines
ax.LineWidth = 1.2;                % Thicker axis lines
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX format

% Add a legend:
legend({'Spatial debris density', 'Max density', 'Departure orbit', 'Arrival orbit'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
hold off

% Initial state in ECI:
r_i = [par.rho_i; 0; 0];                 % Initial position vector (ECI) [km]
v_i = [0; sqrt(par.mu/par.rho_i); 0];    % Initial velocity vector (ECI) [km/s]
S_i = [r_i; v_i];                        % Initial state (ECI) [km, km/s]

% Rotation matrix between initial and final inclination:
T_if = [1, 0, 0; 0, cos(par.Delta_i), -sin(par.Delta_i); 0, sin(par.Delta_i), cos(par.Delta_i)];

% Final state in ECI:
r_f = T_if*[par.rho_f; 0; 0];               % Final position (ECI) [km]
v_f = T_if*[0; sqrt(par.mu/par.rho_f); 0];  % Final velocity (ECI) [km/s]
S_f = [r_f; v_f];                           % Final state (ECI) [km, km/s]

%% Ex 3.2 

% Set dimensionless parameters:
para.hi = 800/par.DU;                          % Altitude of departure orbit [-]
para.i_i = 0;                                  % Initial inclination [rad]
para.hf = 1000/par.DU;                         % Altitude of arrival orbit [-]
para.Delta_i = deg2rad(0.75);                  % Inclination change [rad]
para.i_f = deg2rad(par.i_i + par.Delta_i);     % Final inclination [rad]
para.Re = 6378.1366/par.DU;                    % Earth radius [-]
para.rho_i = (par.hi + par.Re)/par.DU;         % Initial distance from Earth [-]
para.rho_f = (par.hf + par.Re)/par.DU;         % Final distance from Earth [-]
para.mu = 1;                                   % Earth gravitational parameter [-]
para.rho_0 = (750 + par.Re)/par.DU;            % Reference radius for debris flux [-] 
para.k1 = par.k1;                              % Debris spatial density constant 1 [-]
para.k2 = par.k2;                              % Debris spatial density constant 2 [-]
para.m0 = 1000/par.MU;                         % Initial mass [-]
para.T_max = 3/par.FU;                         % Maximum thrust [-]
para.Isp = 3120/par.TU;                        % Specific impulse [-]
para.g0 = par.g0/par.AU;                       % Sea level gravity acceleration [-]

% Dimensionless initial and final state:
Sa_i = [S_i(1:3)/par.DU; S_i(4:6)/par.VU]; % Dimensionless initial state [-]
Sa_f = [S_f(1:3)/par.DU; S_f(4:6)/par.VU]; % Dimensionless final state [-]

% Memorise the variables:
para.Sa_i = Sa_i;
para.Sa_f = Sa_f;

%% Ex 3.4:
% Choose the approach value:
%  1) If approach = 1, the full while cycle is run, so that, starting from 
%     casual solution, the converging solution is chosen. This process is 
%     followed by a solution refinement.
%  2) If approach = 2, there is only solution refinement to speed up the
%     checking process.
approach = 2;

% Define the base part of the initial solution (initial state and mass):
x0_base = [para.Sa_i; para.m0];     % Initial state and mass (dimensionless)
tf0 = 20*pi;                        % Initial time guess (dimensionless)

% Define the options for the propagation:
options_ode = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3), 1e-13*ones(1,1), 1e-13*ones(1,7)]);

% Run full cycle (approach = 1), or pre-loaded solution (approach = 2)
switch approach
    
    case 1 % Run the full cyle
    % Define the options for the optimization:
    options_fsolve = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
        'Display', 'iter', 'MaxFunctionEvaluations', 1e7, ...
        'MaxIterations', 1e6, 'UseParallel', true, ...
        'FiniteDifferenceType','central', ...
        'TolX', 1e-10, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    
    % Useful variable for the cycle:
    attempt = 0;                 % Number of attempt
    Nmax = 100;                  % Maximum number of attempt (to avoid getting stuck)
    success = false;             % Success indicator
    x0 = zeros(8,1);             % Initial guess
    
    % Set always the same random value generation:
    rng default
    
    while ~success && attempt<Nmax
        % Update the attempt number
        attempt = attempt + 1;

        % Construct the guess 
        lambda0_guess(1:6, 1) = -250 + 500*rand(6,1);       % Casual number in the range -250<x<250 are considered.
        lambdam_guess = 250*rand(1,1);                      % Only positive lamda are considered for the one related to mass, as known from theory
        tf0_guess = tf0 -pi + 2*pi*rand(1,1);               % A range of +-pi is considered with respect to the initial guess
        x0 = [lambda0_guess; lambdam_guess; tf0_guess];     % Final guess to be optimised
    
        % Run fsolve:
        [x0_sol, fval_sol, exit_flag] = fsolve(@(x)zerofun([x0_base; x], para, options_ode), x0, options_fsolve);
    
        % Check the optimization results: 
        if exit_flag > 0
            % Stop the while cycle:
            success = true;

        else % In case the attempt fails, print the failure and restart the cycle
            disp(['Attempt number ', num2str(attempt), ' failed. Restart the algorithm with a new guess.']);

        end
    end

    case 2 % Only the results of the previous optimization are saved, to avoid the code run:
        % Lambda values:
        x0_sol = [-214.981220880936; -10.3658727373442; 0.885567773487489; -10.3929207361806; -214.610451351762; -112.94535605739; 2.5964491697761; 64.4801062438446];

        % Zero-finding problem values:
        fval_sol =[-5.67323965583455e-13; -8.21923779037448e-12; 6.87386345760155e-14; 6.95727503385868e-12; -8.42992342597881e-13;  -3.48231339597493e-12; -1.38570188309983e-12 ; -2.11812529321564e-11];
end 

% Propagate the solution found both in state and costate dynamics:
[tt, xx] = ode78(@(t,x) H_dynamics(t, x, para), [0 x0_sol(end)], [x0_base; x0_sol(1:7)], options_ode);

% Initialize the variables:
H = zeros(size(xx,1), 1);              % Hamiltonian function
alpha_u = zeros(size(xx,1), 3);        % Primer vector in ECI frame
alpha_u_NTW = zeros(size(xx,1), 3);    % Primer vector in NTW frame

% Compute the time dependency of the Hamiltonian function and primer vector
% in NTW frame::
for i = 1:size(xx, 1) % Cycle all over the time instants
    
    % Compute the hamiltonian and extract the alpha values:
    [H(i), ~, alpha_u(i, :)] = Hamiltonian(xx(i, :)', para);

    % Transform the control direction in to the NTW frame:
    alpha_u_NTW (i, :) = NTW(alpha_u(i, :), xx(i, 1:3), xx(i, 4:6));

end

% 3D Orbit Plot (with Earth, real dimension)
figure;
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]); % Blue color
grid on;
hold on;
Earth3D(para.Re)
plot3(para.Sa_i(1), para.Sa_i(2), para.Sa_i(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial point'); % Initial point
plot3(para.Sa_f(1), para.Sa_f(2), para.Sa_f(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final point');   % Final point
title(sprintf('Non scaled optimised trajectory in @Earth J2000 for T = %.2f N', 3), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('x [DU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [DU]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [DU]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Trajectory', 'Initial point', 'Final point'}, 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'Box', 'on', 'GridAlpha', 0.3);
set(gcf, 'Color', 'w');
view(45, 30); 
hold off;

% 3D Orbit Plot (without Earth, unreal dimension, but focus on the trajectory)
figure;
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]); % Blue color
grid on;
hold on;
plot3(para.Sa_i(1), para.Sa_i(2), para.Sa_i(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial point'); % Initial point
plot3(para.Sa_f(1), para.Sa_f(2), para.Sa_f(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final point');   % Final point
title(sprintf('Non scaled optimised trajectory in @Earth J2000 for T = %.2f N', 3), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('x [DU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [DU]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [DU]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Trajectory', 'Initial point', 'Final point'}, 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'Box', 'on', 'GridAlpha', 0.3);
set(gcf, 'Color', 'w');
view(45, 30); 
hold off;


% Plot Hamiltonian:
figure;
plot(tt, H, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]); % Blue color with thicker line
grid on;
hold on;
yline(0, '--k', 'LineWidth', 1);
title(sprintf('Hamiltonian Evolution Over Time for T=%.2f N', 3), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('\textbf{Time}[TU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{Hamiltonian Value} $\mathcal{H}(t)$ [-]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Hamiltonian $\mathcal{H}$', 'Null value'}, ...Z
       'FontSize', 14, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');
ylim([-1e-8 1e-8])
box on;
hold off;


% Plot all alpha_u components expressed in NTW reference system:
figure;
plot(tt, alpha_u_NTW(:, 1), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);   % Normal component
hold on;
plot(tt, alpha_u_NTW(:, 2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);        % Transverse component
plot(tt, alpha_u_NTW(:, 3), 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);   % Along-track component
grid on;
title(sprintf('\mathbf{\hat{\alpha}_{u}^*}$ components in NTW frame time evolution for T = %.2f N', 3), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('\textbf{Time} [TU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{$\hat{\alpha}_{u}^*$} [-]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'$\hat{\alpha}_{u, N}^*$ (Normal)', '$\hat{\alpha}_{u, T}^*$ (Transverse)', '$\hat{\alpha}_{u, W}^*$ (Along-track)'}, ...
       'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');
title(sprintf('$\\mathbf{\\hat{\\alpha}_{u}^*}$ Components Evolution Over Time in NTW frame for T = %.2f N', para.T_max * par.FU), ...
    'FontSize', 10, 'Interpreter', 'latex');
hold off;


% Allocate the final results:
solutions.T3.tf = tt(end)*par.TU/60;                 % Final time [min]
solutions.T3.mf = xx(end, 7)*par.MU;                 % Final mass [kg]
solutions.T3.lambdas = x0_sol(1:7);                  % Solving co-state
solutions.T3.err_pos = norm(fval_sol(1:3))*par.DU;   % Error on the position [km]
solutions.T3.err_vel = norm(fval_sol(4:6))*par.VU;   % Error on the velocity [km/s]

% Display the results:
disp('The vector of initial lambdas are');
fprintf('[%.4f] \n', solutions.T3.lambdas);
fprintf('The final time is t_f = %.4f [min] \n', solutions.T3.tf);
fprintf('The final mass is m_f = %.4f [kg] \n', solutions.T3.mf)
fprintf('The final condition error on the position is: %.4e [km] \n', solutions.T3.err_pos);
fprintf('The final condition error on the velocity is: %.4e [km/s] \n', solutions.T3.err_vel);



%% 3.5 Numerical continuation:

% Set the following parameters:
T_fin = 2.860 / par.FU;                      % Adimensionalised final thrust
T_vect = linspace(para.T_max, T_fin, 10);    % Vector for the numerical continuation
x0_old = x0_sol;                             % Initialise the solution from the one computed at the previosu point

% Set the options for the numerical continuation fsolve:
options_fsolve_nc =  optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg', 'Display', 'none', ...
                                  'FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

for i = 2:length(T_vect) % Iterate over the changing thrust
    % Reassign the value of the thrust:
    para.T_max = T_vect(i);

    % Run optimization for decreasing value of thrust, and each time 
    % optimize starting from the value of the previous iterate:
    [x0_new, fval_new, exit_flag] = fsolve(@(x)zerofun([x0_base; x], para, options_ode), x0_old, options_fsolve_nc);
 
    % Transform the new variable into the new to restart the cycle:
    if exit_flag > 0 % If the algorithm reaches convergence
        x0_old = x0_new; % Re-allocate the initial guess

    else % The algorithm did not reach convergence
        warning('The numerical continuation for T = %.4f [N] did not reach convergence', T_vect(i)*par.FU);
    end
end

% Propagate the solution found after the last optimization:
[tt_last, xx_last] = ode78(@(t,x) H_dynamics(t, x, para), [0 x0_new(end)], [x0_base; x0_new(1:7)], options_ode);

% Initialise the variables
H_last = zeros(size(xx_last,1), 1);            % Hamiltonian function for the desired level of thrust
alpha_u_last = zeros(size(xx_last,1), 3);      % Primer vector components in cartesian state in Earth J200
alpha_u_NTW_last = zeros(size(xx_last,1), 3);  % Primer vector components in NTW frame

% Compute the time dependency of the Hamiltonian function and primer vector
% in NTW frame:
for i = 1:size(xx_last, 1)
    
    % Compute the Hamiltonian and extract the alpha values:
    [H_last(i), ~, alpha_u_last(i, :)] = Hamiltonian(xx_last(i, :)', para);

    % Transform the control direction into the NTW frame:
    alpha_u_NTW_last(i, :) = NTW(alpha_u_last(i, :), xx_last(i, 1:3), xx_last(i, 4:6));
end

% 3D Orbit Plot from Last Optimization (with Earth and correct scaling):
figure;
plot3(xx_last(:, 1), xx_last(:, 2), xx_last(:, 3), 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]); % Purple color
grid on;
hold on;
Earth3D(para.Re)
plot3(para.Sa_i(1), para.Sa_i(2), para.Sa_i(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial Point');
plot3(para.Sa_f(1), para.Sa_f(2), para.Sa_f(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final Point');
title(sprintf('Non scaled optimised trajectory in @Earth J2000 for T = %2.f N', T_fin * par.FU), 'FontSize', 14, 'Interpreter', 'latex');
xlabel('x [DU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [DU]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [DU]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Trajectory', 'Initial Point', 'Final Point'}, 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'Box', 'on', 'GridAlpha', 0.3);
set(gcf, 'Color', 'w'); 
view(45, 30); 
hold off;

% Plot the 3D orbit (without Earth, to highligth the orbits themselves)
figure;
plot3(xx_last(:, 1), xx_last(:, 2), xx_last(:, 3)*40, 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]); % Purple color
grid on;
hold on;
plot3(para.Sa_i(1), para.Sa_i(2), para.Sa_i(3)*40, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial Point');
plot3(para.Sa_f(1), para.Sa_f(2), para.Sa_f(3)*40, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final Point');
title(sprintf('(Scaled) 3D Trajectory in ECI for T = %4.f N', T_fin * par.FU), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('x [DU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [DU]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [DU]', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Trajectory', 'Initial Point', 'Final Point'}, 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'Box', 'on', 'GridAlpha', 0.3);
set(gcf, 'Color', 'w'); 
view(45, 30); 
hold off;

% Plot Hamiltonian from Last Optimization:
figure;
plot(tt_last, H_last, 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]); % Purple color
grid on;
hold on;
yline(0, '--k', 'LineWidth', 1); % Reference line at y = 0
title(sprintf('Hamiltonian evolution over time for T = %.2f N', T_fin * par.FU), 'FontSize', 14, 'Interpreter', 'latex');
xlabel('\textbf{Time} [TU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{Hamiltonian Value} $\mathcal{H}(t) [-]$', 'FontSize', 12, 'Interpreter', 'latex');
legend({'Hamiltonian $\mathcal{H}$', 'Null value'}, 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex');
ylim([-10e-9 10e-9])
set(gcf, 'Color', 'w'); % White background
hold off;

% Plot all alpha_u components expressed in NTW reference system:
figure;
plot(tt_last, alpha_u_NTW_last(:, 1), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);   % Normal component
hold on;
plot(tt_last, alpha_u_NTW_last(:, 2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);        % Transverse component
plot(tt_last, alpha_u_NTW_last(:, 3), 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);   % Along-track component
hold off;
grid on;

% Fixing title syntax
title(sprintf('$\\mathbf{\\hat{\\alpha}_{u}^*}$ Components Evolution Over Time in NTW frame for T = %.2f N', T_fin * par.FU), ...
    'FontSize', 10, 'Interpreter', 'latex');

% Fixing xlabel and ylabel syntax
xlabel('\textbf{Time} [TU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel(sprintf('$\\mathbf{\\hat{\\alpha}_{u}^*}$ [-]'), 'FontSize', 12, 'Interpreter', 'latex');

% Fixing legend syntax
legend({'$\hat{\alpha}_{u, N}^*$ (Normal)', '$\hat{\alpha}_{u, T}^*$ (Transverse)', '$\hat{\alpha}_{u, W}^*$ (Along-track)'}, ...
       'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
% Setting properties for the plot
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

% Plot the 3D orbit comparison (without Earth, to highligth the orbits themselves)
figure;
plot3(xx_last(:, 1), xx_last(:, 2), xx_last(:, 3)*40, 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', sprintf('Optimised trajectory with T = %.2f N', par.T_max)); % Purple color
grid on;
hold on;
plot3(xx(:, 1), xx(:, 2), xx(:, 3)*40, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName', sprintf('Optimised trajectory with T = %.2f N', T_fin*par.FU)); % Blue color
plot3(para.Sa_i(1), para.Sa_i(2), para.Sa_i(3)*40, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial Point');
plot3(para.Sa_f(1), para.Sa_f(2), para.Sa_f(3)*40, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final Point');
title(sprintf('Comparison between the (scaled) 3D Trajectories in @Earth J2000'), 'FontSize', 10, 'Interpreter', 'latex');
xlabel('x [DU]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [DU]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [DU]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'FontSize', 10, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 10, 'Box', 'on', 'GridAlpha', 0.3);
set(gcf, 'Color', 'w'); 
view(45, 30); 
hold off;

% Plot the comparison between the spatial debri densities
figure; hold on; grid on;
plot(tt_last*par.TU/60, q(vecnorm(xx_last(:, 1:3), 2, 2) * par.DU), 'LineWidth', 2, ...
    'Color', [1, 0, 0], 'DisplayName', sprintf('Optimised trajectory with T = %.2f N', T_fin*par.FU));
plot(tt*par.TU/60, q(vecnorm(xx(:, 1:3), 2, 2) * par.DU), 'LineWidth', 2, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', sprintf('Optimised trajectory with T = %.2f N', par.T_max));
legend('show', 'Interpreter' , 'latex', 'FontSize', 14, 'location', 'best');
title('Spatial debris density along the optimised trajectories', 'Interpreter', 'latex');
ylabel('Spatial debris density [$DU^{-3}$]', 'Interpreter','latex', 'FontSize', 18);
xlabel('Time along the trajectory [min]', 'Interpreter','latex', 'FontSize', 18);
plot(tt_last(1)*par.TU/60, q(norm(xx_last(1, 1:3)) * par.DU), 'DisplayName', 'Initial Orbit', ...
    'Marker','o', 'Color', [0 1 0], 'MarkerSize', 8, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'k');
plot(tt_last(end)*par.TU/60, q(norm(xx_last(end, 1:3)) * par.DU), 'DisplayName', ...
    'Final Orbit','Marker','o', 'Color', [0 0 1], 'MarkerSize', 8, 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', 'k');
hold off;


% Allocate final results:
solutions.T2.tf = tt_last(end)*par.TU/60;            % Final time [min]
solutions.T2.mf = xx_last(end, 7)*par.MU;            % Final mass [kg]
solutions.T2.lambdas = x0_new(1:7);                  % Solving co-states for T_fin
solutions.T2.err_pos = norm(fval_new(1:3))*par.DU;   % Error on position in ECI [km]
solutions.T2.err_vel = norm(fval_new(4:6))*par.VU;   % Error on velocity in ECI [km/s]

% Display the results:
fprintf('The vector of initial lambdas is \n');
fprintf('[%.4f]\n', solutions.T2.lambdas);
fprintf('The final time is t_f = %.4f min \n', solutions.T2.tf);
fprintf('The final mass is m_f = %.4f kg \n', solutions.T2.mf)
fprintf('The final condition error on the position is: %.4e [km] \n', solutions.T2.err_pos);
fprintf('The final condition error on the velocity is: %.4e [km/s] \n', solutions.T2.err_vel);

%% Functions:
%--------------------------------------------------------------------------
function dx = H_dynamics(~, x, para)
% DESCRIPTION:
% Computes the state and costate dynamics based on the Hamiltonian 
% dynamics in an optimal control problem.
%
% PROTOTYPE:
% dx = H_dynamics(~, x, para)
%
% INPUT:
% ~       Time variable (placeholder, not used explicitly)          [-]
% x [14x1] State vector containing: 
%             - Position vector rr (x(1:3))                         [-]
%             - Velocity vector vv (x(4:6))                         [-]
%             - Mass m (x(7))                                       [-]
%             - Costate of position llr (x(8:10))                   [-]
%             - Costate of velocity llv (x(11:13))                  [-]
%             - Lagrangian multiplier (x(14))                       [-]
% para    Struct containing problem parameters: 
%             - mu       Earth gravitational parameter              [-]
%             - T_max    Maximum thrust                             [-]
%             - Isp      Specific impulse                           [-]
%             - g0       Gravity acceleration at sea level          [-]
%             - k1       Lagrangian function parameter 1            [-]
%             - k2       Lagrangian function parameter 2            [-]
%             - rho_0    Initial distance from Earth                [-]
%
% OUTPUT:
% dx[14x1] Derivative of the state vector, including dynamics for:
%             - Position, velocity, mass, costate of position and velocity.
%
    % Some useful parameters:
    mu = para.mu;               % Earth gravitational parameter [-]
    T_max = para.T_max;         % Maximum thrust [-]
    Isp = para.Isp;             % Specific impulse [-]
    g0 = para.g0;               % Sea level gravity acceleration [-]
    k1 = para.k1;               % Lagrangian function parameter 1
    k2 = para.k2;               % Lagrangian function parameter 2
    rho0_DU = para.rho_0;       % Initial distance from Earth [-]


    % Rename conveniently the variables of final state:
    rr = x(1:3);                            % Position vector
    r = norm(rr);                           % Norm of the position vector
    vv = x(4:6);                            % Velocity vector
    m = x(7);                               % Mass
    llr = x(8:10);                          % Costate associated with the position
    llv = x(11:13);                         % Costate associated to the velocity
    lv = norm(llv);                         % Norm of the costate associated with the velocity

    % Basing on the solution of PMP, you get that the optmimal
    % thrust-pointing direction and magnitude are:
    alpha_u = -llv/lv;    % Optimal thrust-pointing direction
    u = 1;                % Control magnitude

    % Allocate the space for the derivative of the state:
    dx = zeros(14,1);

    % Impose the state dynamics
    dx(1:7) = [vv;                                 % Velocity state dynamics
               -(mu/r^3)*rr + u*T_max*alpha_u/m;   % Acceleration state dynamics
               -u*T_max/(Isp*g0)];                 % Mass state dynamics
    
    % Impose the costate dynamics
    dq_dr = 2*k1*(r-rho0_DU)/((k2+(r-rho0_DU)^2)^2*r).*rr;        % Derivative of lagrangian function with respect to the position
    dx(8:14) =  [-3*mu/r^5*dot(rr,llv)*rr + mu/r^3*llv + dq_dr;   % Position costate dynamics
                -llr;                                             % Velocity costate dynamics
                -u*lv*T_max/m^2];                                 % Mass costate dynamics

end

%--------------------------------------------------------------------------
function [H, f, alpha_u] = Hamiltonian(x, para)
% DESCRIPTION:
% Computes the Hamiltonian (H), the state dynamics (f), and the 
% optimal thrust-pointing direction (alpha_u) based on the given 
% state and problem parameters.
%
% PROTOTYPE:
% [H, f, alpha_u] = Hamiltonian(x, para)
%
% INPUT:
% x[14x1]      State vector containing:
%                - Position vector rr (x(1:3))                     [-]
%                - Velocity vector vv (x(4:6))                     [-]
%                - Mass m (x(7))                                   [-]
%                - Costate of position llr (x(8:10))               [-]
%                - Costate of velocity llv (x(11:13))              [-]
%                - Costate of mass lm (x(14))                      [-]
% para         Struct containing problem parameters:
%                - mu      Earth gravitational parameter            [-]
%                - rho_0   Initial distance from Earth              [-]
%                - T_max   Maximum thrust                           [-]
%                - Isp     Specific impulse                         [-]
%                - g0      Gravity acceleration at sea level        [-]
%                - k1, k2  Spatial density constants                [-]
%
% OUTPUT:
% H[1x1]       Hamiltonian value                                   [-]
% f[7x1]       State dynamics                                      [-]
% alpha_u[3x1] Optimal thrust-pointing direction                   [-]
%

    % Extract from para the useful adimensional parameters:
    mu = para.mu;               % Earth gravitational parameter [-]
    rho0_DU = para.rho_0;       % Initial distance from Earth [-]
    T_max = para.T_max;         % Maximum thrust [-]
    Isp = para.Isp;             % Specific impulse [-]
    g0 = para.g0;               % Sea level gravity acceleration [-]
    k1 = para.k1;               % Debris spatial density constant 1 [-]
    k2 = para.k2;               % Debris spatial density constant 2 [-]

    % Rename conveniently the variables of final state:
    rr = x(1:3);                % Position vector
    r = norm(rr);               % Norm of the position vector            
    vv = x(4:6);                % Velocity vector
    m = x(7);                   % Mass
    llr = x(8:10);              % Costate associated with the position
    llv = x(11:13);             % Costate associated with the velocity
    lv = norm(llv);             % Norm of the costate associated with the velocity vector
    lm = x(14);                 % Costate associated with the mass
    ll = [llr; llv; lm];        % (Full) Costate vector

    % Basing on the solution of PMP, you get that the optmimal
    % thrust-pointing direction and magnitude are:
    alpha_u = -llv/lv;
    u = 1;
    
    % Define the Hamiltonian function in the final state:
    l = k1/(k2 + (r-rho0_DU)^2);                 % Objective function
    f = [vv;                                     % Velocity ynamics
         -mu/r^3*rr + u*T_max*alpha_u/m;         % Acceleration dynamics
         -u*T_max/(Isp*g0)];                     % Mass dynamics
    H = l + dot(ll, f);                          % Hamiltonian function
end

%--------------------------------------------------------------------------
function fun = zerofun(AS, para, options_STM)
% DESCRIPTION:
% Computes the final condition residuals for a given state and 
% propagates the dynamics to satisfy the required constraints.
%
% PROTOTYPE:
% fun = zerofun(AS, para, options_STM)
%
% INPUT:
% AS[nx1]       Augmented state vector, including:
%                 - Initial state S0 (AS(1:end-1))                 [-]
%                 - Final time tf (AS(end))                        [-]
% para          Struct containing problem parameters:
%                 - Sa_f Final desired state vector                [-]
% options_STM   ODE solver options for state transition matrix 
%               (default: relative tolerance = 1e-12)              [-]
%
% OUTPUT:
% fun[8x1]      Residuals of the final conditions                  [-]
%

    % Allocate the variables:
    S0 = AS(1:(end-1),1);  % Initial state 
    t0 = 0;                % Initial time
    tf = AS(end,1);        % Final time
    
    % Extract from para the useful adimensional parameters:
    Sa_f = para.Sa_f;

    if nargin < 3 % If not specified, assign automatically a reasonable value of the tolerance:
        options_STM = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-13*ones(1,3), 1e-13*ones(1,1), 1e-13*ones(1,7)]);
    end
    
    % Propagate up to the final conditions:
    [~, xx] = ode78(@(t,x) H_dynamics(t, x, para), [t0 tf], S0, options_STM);

    % Extract (final) state vector:
    xf = xx(end, 1:14)';

    % Calculate the Hamiltonian and the RHS of the dynamics of the final 
    % state:
    [Hf, ~, ~] = Hamiltonian(xf, para);
 
    % Impose the condition to be minimised:
    fun(1:3) = Sa_f(1:3) - xf(1:3);           % Impose final position
    fun(4:6) = Sa_f(4:6) - xf(4:6);           % Impose final velocity
    fun(7) = xf(14);                          % Condition on mass costate
    fun(8) = Hf;                              % Trasversality condition

end

%--------------------------------------------------------------------------
function alpha_NTW = NTW(alpha, r, v)
% DESCRIPTION:
% Transforms a control vector into the NTW reference frame.
%
% PROTOTYPE:
% alpha_NTW = NTW(alpha, r, v)
%
% INPUT:
% alpha[1x3]    Control vector in the inertial frame            [-]
% r[1x3]        Position vector                                 [-]
% v[1x3]        Velocity vector                                 [-]
%
% OUTPUT:
% alpha_NTW[1x3] Control vector in the NTW reference frame      [-]
%

    % Define the NTW reference frame
    T_axis = v / norm(v);                        % Tangential axis
    N_axis = cross(r, T_axis);                   % Normal axis (pre-normalization)
    N_axis = N_axis / norm(N_axis);              % Normalize the normal axis
    W_axis = cross(N_axis, T_axis);              % Out-of-plane axis

    % Error checking for NTW orthogonality
    if abs(dot(N_axis, T_axis)) > 1e-12 || ...   % Check orthogonality between N and T
       abs(dot(N_axis, W_axis)) > 1e-12          % Check orthogonality between N and W
        disp('Error in the NTW reference system definition: Axes are not orthogonal.')
        return;
    elseif abs(norm(cross(N_axis, cross(r, v) / norm(cross(r, v))))) > 1e-12
        % Check if N-axis is parallel to the angular momentum
        disp('N-axis is not parallel to the angular momentum.')
    end 

    % Transformation matrix for NTW frame
    T = [N_axis;                                 % First row: Normal axis
         T_axis;                                 % Second row: Tangential axis
         W_axis];                                % Third row: Out-of-plane axis

    % Transform the control vector to the NTW reference system
    alpha_NTW = T * alpha';                      % Compute control vector in NTW frame

end

%--------------------------------------------------------------------------
function Earth3D(Rt)
% DESCRIPTION
% This function plots the Earth in the center of the 3D plot with the
% radius Rt specified
%
% PROTOTYPE
% Earth3D(Rt)
%
% INPUT:
%  Rt [1]     Radius to be assigned to the Earth plot [km]
%             If not specified, Rt is automatically set to 6378.1363 km.
%
% OUTPUT:
%  Plot
    
    if nargin == 0 % If no input are passed
        Rt = 6378.1363; % Default Earth's radius in km
    end
    
    % Load topographic data
    load('topo.mat', 'topo', 'topomap1');
    
    % Generate the spherical surface
    [x, y, z] = sphere(100);
    
    % Scale the sphere by Earth's radius
    x = x * Rt;
    y = y * Rt;
    z = z * Rt;
    
    % Apply topography as texture
    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = 0.5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor = 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = flip(topo, 1); % Flip the texture to align correctly
    h = surface(x, y, z , props);
    
    % Turn off handle visibility to exclude from legend
    set(h, 'HandleVisibility', 'off');
    
    % Add lighting
    light('Position', [1 1 1], 'Style', 'infinite');
    light('Position', [-1.5 0.5 -0.5], 'Style', 'local', 'Color', [0.6 0.2 0.2]);
    
    axis equal; % Ensure equal scaling
    hold on;    % Impose hold on to conitnue plotting also outside the function

end
