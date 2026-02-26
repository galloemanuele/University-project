% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #1, Exercise 2
% Author: Emanuele Gallo

% Prepare the workspace:
clear 
close all
clc

% Impose long format to better visualize results:
format long g

% Load the kernels:
cspice_furnsh('ex02.tm')

% Extract the useful data:
par.E_Radii = cspice_bodvrd('Earth', 'RADII',3);             % radius of Earth
par.E_GM = cspice_bodvrd('Earth', 'GM', 1);                  % for Earth
par.M_Radii = cspice_bodvrd('Moon', 'RADII',3);              % radius of the Moon
par.M_GM = cspice_bodvrd('Moon', 'GM', 1);                   % for moon
par.mu = par.M_GM/(par.E_GM+par.M_GM);                       % gravitational parameter

% Constants to solve the PCRTBP:
par.ms = 3.28900541e5;      % [-]: Scaled mass of the Sun
par.rho = 3.88811143e2;     % [-]: Scaled Sun-(Earth+Moon) distance
par.ws = -9.25195985e-1;    % [-]: Scaled angular velocity of the Sun
par.wem = 2.66186135e-1;    % [s^(-1)]: Earth-Moon angular velocity 
par.lem = 3.84405e8;        % [m]: Earth-Moon distance
par.hi = 167;               % [km]: Altitude of departure orbit
par.hf = 100;               % [km]: Altitude of arrival orbit
par.DU = 3.84405000e5;      % [km]: Distance Unit
par.TU = 4.34256461;        % [days]: Time Unit
par.VU = 1.02454018;        % [km/s]: Velocity Unit
par.T_orbit = 2*pi;         % [-]: Orbital period (adimensional)

% Parking orbit radius initial and final:
par.ri = (par.E_Radii(1)+par.hi)/par.DU;
par.rf = (par.M_Radii(1)+par.hf)/par.DU;


%% Exercise 2.1

% Variables to construct the initial guess:
alpha0 = 0.2*pi;            % angle defined on the Earth circular parking orbit
beta0 = 1.41;               % initial-to-circular velocity ratio
delta0 = 4;                 % transfer duration
ti0 = 2;                    % initial time

% Calculate the first guess solution using a custom function:
AS0 = first_guess(alpha0, beta0, delta0, ti0, par);    % Generates first guess for trajectory and time
tf0 = AS0(end);                                        % Extract the final time from the guess solution

% Calculate the first guess solution:
fprintf('-------------------------First guess solution------------------------------------------\n');
disp('The optimised initial state is:')
fprintf('x = %.6f [DU] \n', AS0(1));
fprintf('y = %.6f [DU] \n', AS0(2));
fprintf('vx = %.6f [VU] \n', AS0(3));
fprintf('vy = %.6f [VU] \n', AS0(4));
fprintf('The initial time is: %.3f [TU] \n', AS0(5));
fprintf('The final time is: %.3f [TU] \n', AS0(6));
fprintf('The corresponding time of flight is: %.3f [TU]\n', AS0(6)-AS0(5));  % Display the time of flight
fprintf('Quantities expressed @EMB Earth-Moon rotating frame. \n \n \n ');

% Propagate in the rotating frame, thus calculating the trajectory of the 
% transfer and the related times:
[~, PHIf, tf, xx_trajectory0, tt_trajectory0]  = propagation(ti0, AS0(1:4), tf0, par, 'xyPBRFBP_STM_ROT');

% Convert in the Earth-centered reference frame:
xx_EC = ROT2ECI(xx_trajectory0(:, 1:4), tt_trajectory0, par); % The first guess solution
xx_M_EC = ROT2ECI([(1-par.mu)*ones(size(tt_trajectory0, 1), 1), zeros(size(tt_trajectory0,1), 3)], tt_trajectory0 , par);

% Plot for trajectory in rotating reference frame (first subplot)
figure
plot(xx_trajectory0(:,1), xx_trajectory0(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
hold on
grid on
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8) % Moon
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8) % Earth
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
axis equal
title('Trajectory in the Earth-Moon rotating frame centred @EMB', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Trajectory', 'Moon', 'Earth'}, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14)
hold off

% Plot for trajectory in Earth-centered reference frame (second subplot)
figure
plot(xx_EC(:,1), xx_EC(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
hold on
plot(xx_M_EC(:,1), xx_M_EC(:,2),'r--','LineWidth', 1.5)
plot(xx_M_EC(end,1), xx_M_EC(end,2), 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8) % Moon
grid on
plot(0, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8) % Earth
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
axis equal
title('Trajectory in the Earth-centered Inertial reference frame', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Trajectory', 'Moon trajectory', 'Moon position at final instant', 'Earth'}, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14)
hold off


%% Exercise 2.2: Simple shooting

% The lower and upper boundaries are set from theoretical and physical 
% considerations:
v_margin = 1.5*sqrt(par.E_GM/(par.ri*par.DU));                                        % Margins on velocity
lb_s = [-par.ri, -par.ri, -v_margin, -v_margin, 0, 0];                                % Lower bounds
ub_s = [par.ri, par.ri, v_margin, v_margin, -2*pi/par.ws, -2*pi/par.ws + 28/par.TU];  % Upper bounds

% Set the linear inequality constraint (ti-tf<0):
A = [zeros(1, 4), 1, -1]; % Linear inequality constraint matrix
b = 0;                    % Linear inequality constraint known vector

% Perform the optimization without providing the analytical form of the 
% gradient:
options_NG = optimoptions('fmincon', 'Algorithm', 'active-set', ...                          % Choose active set algorithm to solve
                   'Display', 'none', ...                                                    % Show detailedly what the algorithm does at each iterate
                   'MaxIter', 100, ...                                                       % Set maximum number of iterations
                   'TolFun', 1e-12, ...                                                      % Absolute tolerance on the function value
                   'ConstraintTolerance', 1e-10, ...
                   'GradObj', 'off', ...
                   'GradCon', 'off');  
[AS0_opt_ng_s, f_ng_s, exflag_ng_s, output_ng_s] = fmincon(@(ASi)ObjFun(ASi, par, 'simple'), AS0, A, b, [], [], ...
                                           lb_s, ub_s, @(ASi)constraints(ASi, par, 'simple'), options_NG); 

% Display the solution:
fprintf('-------------------------Simple shooting method without gradient-----------------------------------\n');
disp('The optimised initial state is:')
fprintf('x = %.6f [DU] \n', AS0_opt_ng_s(1));
fprintf('y = %.6f [DU] \n', AS0_opt_ng_s(2));
fprintf('vx = %.6f [VU] \n', AS0_opt_ng_s(3));
fprintf('vy = %.6f [VU] \n', AS0_opt_ng_s(4));
fprintf('The initial time is: %.3f [TU] \n', AS0_opt_ng_s(5));
fprintf('The final time is: %.3f [TU] \n', AS0_opt_ng_s(6));
fprintf('The corresponding time of flight is: %.3f [TU] \n', AS0_opt_ng_s(6)-AS0_opt_ng_s(5));  % Display the time of flight
fprintf('The optimised value of the objective function (i.e. the cost) is: %.4f [VU]\n', f_ng_s);  % Display the objective function value
fprintf('Quantities expressed @EMB Earth-Moon rotating frame. \n \n \n ')

% Propagation of the found solution:
[Sf_ng_s, ~, ~, xx_ng, tt_ng]  = propagation(AS0_opt_ng_s(5), AS0_opt_ng_s(1:4), AS0_opt_ng_s(6), par, 'xyPBRFBP_STM_ROT');

% Characterise the final orbit:
[orbitType_ng_s, captureType_ng_s] = characterise_orbits(Sf_ng_s, par);
disp('The achieved orbit can be characterised as follows:') % Display the orbit and capture types
fprintf(' - %s \n', orbitType_ng_s); 
fprintf(' - %s \n \n', captureType_ng_s)  

% Plot the optimized trajectory without gradient in the rotating frame
figure
plot(xx_ng(:,1), xx_ng(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
hold on
grid on

% Plot the circular orbit around Earth (radius par.ri)
theta = linspace(0, 2*pi, 100); % Parametric angle for a full circle
earth_orbit_x = -par.mu + par.ri * cos(theta); % x coordinates of Earth orbit
earth_orbit_y = par.ri * sin(theta); % y coordinates of Earth orbit
plot(earth_orbit_x, earth_orbit_y, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2) % Earth orbit

% Plot the circular orbit around the Moon (radius par.rf)
moon_orbit_x = 1 - par.mu + par.rf * cos(theta); % x coordinates of Moon orbit
moon_orbit_y = par.rf * sin(theta); % y coordinates of Moon orbit
plot(moon_orbit_x, moon_orbit_y, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2) % Moon orbit

% Plot the Moon and Earth positions
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8) % Moon
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8) % Earth

% Labels and formatting
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
legend({'Trajectory', 'Orbit around the Earth', 'Orbit around the Moon', 'Moon', 'Earth'}, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10)
title('Trajectory optimised without gradient (Earth-Moon rotating frame centred @EMB)', 'FontSize', 10)
axis equal
% Add zoomed subplots
% Near Earth
axes('Position', [0.6, 0.5, 0.23, 0.23]); % Position zoom at bottom-right
plot(xx_ng(:,1), xx_ng(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);  % Plot trajectory near Earth
hold on
plot(earth_orbit_x, earth_orbit_y, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2);  % Earth orbit
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8);  % Earth position
axis equal
xlim([-1.2*par.ri, 1.2*par.ri] - par.mu);  % Adjust x-axis limits
ylim([-1.2*par.ri, 1.2*par.ri]);  % Adjust y-axis limits
grid on;
box on;
set(gca, 'xtick', [], 'ytick', []); % Remove ticks

% Near Moon
axes('Position', [0.2, 0.5, 0.23, 0.23]); % Position zoom at bottom-left
plot(xx_ng(:,1), xx_ng(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);  % Plot trajectory near Moon
hold on
plot(moon_orbit_x, moon_orbit_y, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);  % Moon orbit
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8);  % Moon position
axis equal
xlim([-1.2*par.rf, 1.2*par.rf] + (1 - par.mu));  % Adjust x-axis limits
ylim([-1.2*par.rf, 1.2*par.rf]);  % Adjust y-axis limits
grid on;
box on;
set(gca, 'xtick', [], 'ytick', []);  % Remove ticks
hold off



% Set the options for gradient optimisation:
options_G_s = optimoptions('fmincon', 'Algorithm', 'active-set', ...                 % Choose active set algorithm to solve
               'Display', 'none', ...                                                % Show detailedly what the algorithm does at each iterate
               'MaxIter', 1000, ...                                                  % Set maximum number of iterations
               'FunctionTolerance', 1e-11, ...                                       % Absolute tolerance on the function value
               'ConstraintTolerance', 1e-11, ...
               'FiniteDifferenceType', 'central', ...
               'FiniteDifferenceStep', 1e-7,...
               'MaxFunctionEvaluations', 1000, ...                                   % Set max number of function evaluations
               'SpecifyObjectiveGradient', true, ...
               'SpecifyConstraintGradient', true, ...
               'CheckGradients', false);    
[AS0_opt_g_s, f_g_s, exflag_g_s, output_g_s] = fmincon(@(ASi)ObjFun(ASi, par, 'simple'), AS0, A, b, [], [], ...
                                         [], [], @(ASi)constraints(ASi, par, 'simple'), options_G_s);

% Display the solution:
fprintf('-------------------------Simple shooting method with gradient---------------------------------------\n');
disp('The optimised initial state is:')
fprintf('x = %.6f [DU] \n', AS0_opt_g_s(1));
fprintf('y = %.6f [DU] \n', AS0_opt_g_s(2));
fprintf('vx = %.6f [VU] \n', AS0_opt_g_s(3));
fprintf('vy = %.6f [VU] \n', AS0_opt_g_s(4));
fprintf('The initial time is: %.3f [TU] \n', AS0_opt_g_s(5));
fprintf('The final time is: %.3f [TU] \n', AS0_opt_g_s(6));
fprintf('The corresponding time of flight is: %.3f [TU]\n', AS0_opt_g_s(6)-AS0_opt_g_s(5));  % Display the time of flight
fprintf('The optimised value of the objective function (i.e. the cost) is: %.4f [VU]\n', f_g_s);  % Display the objective function value
fprintf('Quantities expressed @EMB Earth-Moon rotating frame. \n \n \n ')

% Propagation of the found solution:
[Sf_g_s, ~, ~, xx_g_s, tt_g_s]  = propagation(AS0_opt_g_s(5), AS0_opt_g_s(1:4), AS0_opt_g_s(6), par, 'xyPBRFBP_STM_ROT');

% Characterise the final orbit:
[orbitType_g_s, captureType_g_s] = characterise_orbits(Sf_g_s, par);
disp('The achieved orbit can be characterised as follows:') % Display the orbit and capture types
fprintf(' - %s \n', orbitType_g_s); 
fprintf(' - %s \n \n', captureType_g_s);  

% Plot the optimized trajectory with gradient in the rotating frame
figure
plot(xx_g_s(:,1), xx_g_s(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
hold on
grid on

% Plot the circular orbit around Earth (radius par.ri)
theta = linspace(0, 2*pi, 100); % Parametric angle for a full circle
earth_orbit_x = -par.mu + par.ri * cos(theta); % x coordinates of Earth orbit
earth_orbit_y = par.ri * sin(theta); % y coordinates of Earth orbit
plot(earth_orbit_x, earth_orbit_y, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2) % Earth orbit

% Plot the circular orbit around the Moon (radius par.rf)
moon_orbit_x = 1 - par.mu + par.rf * cos(theta); % x coordinates of Moon orbit
moon_orbit_y = par.rf * sin(theta); % y coordinates of Moon orbit
plot(moon_orbit_x, moon_orbit_y, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2) % Moon orbit

% Plot the Moon and Earth positions
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8) % Moon
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8) % Earth

% Labels and formatting
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)
legend({'Trajectory', 'Orbit around the Earth', 'Orbit around the Moon', 'Moon', 'Earth'}, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10)
axis equal
title('Trajectory optimised with gradient (Earth-Moon rotating frame centred @EMB)', 'FontSize', 10)
% Add zoomed subplots for the optimized trajectory with gradient
% Near Earth
axes('Position', [0.6, 0.5, 0.23, 0.23]);  % Set position of zoomed plot near Earth (bottom-right)
plot(xx_g_s(:,1), xx_g_s(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);  % Plot the trajectory
hold on  % Hold the plot to overlay the next elements
plot(earth_orbit_x, earth_orbit_y, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2);  % Plot Earth orbit
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8);  % Plot Earth position
axis equal  % Set equal scaling for x and y axes
xlim([-1.2*par.ri, 1.2*par.ri] - par.mu);  % Limit x-axis range around Earth orbit
ylim([-1.2*par.ri, 1.2*par.ri]);  % Limit y-axis range around Earth orbit
grid on;  % Turn on the grid
box on;  % Add box around the plot
set(gca, 'xtick', [], 'ytick', []);  % Remove x and y ticks

% Near Moon
axes('Position', [0.2, 0.5, 0.23, 0.23]);  % Set position of zoomed plot near Moon (bottom-left)
plot(xx_g_s(:,1), xx_g_s(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);  % Plot the trajectory
hold on  % Hold the plot to overlay the next elements
plot(moon_orbit_x, moon_orbit_y, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);  % Plot Moon orbit
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8);  % Plot Moon position
axis equal  % Set equal scaling for x and y axes
xlim([-1.2*par.rf, 1.2*par.rf] + (1 - par.mu));  % Limit x-axis range around Moon orbit
ylim([-1.2*par.rf, 1.2*par.rf]);  % Limit y-axis range around Moon orbit
grid on;  % Turn on the grid
box on;  % Add box around the plot
set(gca, 'xtick', [], 'ytick', []);  % Remove x and y ticks
hold off  % End the plot hold



%% Exercise 2.3: Multiple shooting

% Number of shootings:
N = 4;

% Corresponding number elements within an interval:
N_el_interval = floor(length(tt_trajectory0)/(N-1));

% Create the proper variables useful to optimize the problem both for the
% initial guess and the boundaries:
t_vect = zeros(length(N), 1);
lb_m = zeros(1, 4*N+2);
ub_m = zeros(1, 4*N+2);
AS0_m = zeros(4*N, 1);
AS0_m (1:4, 1) = AS0(1:4);

for i=1:N

    % Vector of times:
    t_vect(i) = ti0 + (i-1)/(N-1)*(tf0-ti0);

    % Initial guess for the fmincon:
    if i==1 % Initial 
        % Boundaries
        lb_m(1, 4*(i-1) + 1:4*i) = [-4, -4, -inf, -inf];
        ub_m(1, 4*(i-1) + 1:4*i) = [4, 4, inf, inf];
    
    elseif i>1
        % Propagation to get initial condition to the solver:
        [AS0_m_na, ~, ~, ~, ~]  = propagation(t_vect(i-1), AS0_m(4*(i-2)+1:4*(i-1),1), t_vect(i), par, 'xyPBRFBP_STM_ROT');
        AS0_m(4*(i-1)+1:4*i,1) = AS0_m_na;

        % Boundaries:
        lb_m(1, 4*(i-1) + 1:4*i) = [-(AS0_m_na(1) + 4), -(AS0_m_na(2) + 4), -(AS0_m_na(3) + 5*v_margin), -(AS0_m_na(4) + 5*v_margin)];
        ub_m(1, 4*(i-1) + 1:4*i) = [(AS0_m_na(1) + 4), (AS0_m_na(2) + 4), (AS0_m_na(3) + 5*v_margin), (AS0_m_na(4) + 5*v_margin)];

    end 

end

AS0_m = [AS0_m; ti0; tf0];
lb_m(1, end-1:end) = [0, 2/par.TU];
ub_m(1, end-1:end) = [28/par.TU, inf];

% Solve using fmincon
options_G_m = optimoptions('fmincon', 'Algorithm', 'active-set', ...                            % Choose active set algorithm to solve
                       'Display', 'none', ...                                                   % Show detailedly what the algorithm does at each iterate
                       'MaxIter', 5000, ...                                                     % Set maximum number of iterations
                       'FiniteDifferenceType', 'central', ...
                       'FiniteDifferenceStepSize', sqrt(eps), ...
                       'ConstraintTolerance', 1e-9, ...
                       'MaxFunctionEvaluations', 1e10, ...                                      % Set max number of function evaluations
                       'CheckGradients', false, ...
                       'SpecifyObjectiveGradient', true, ...
                       'SpecifyConstraintGradient', true); 

[AS0_opt_g_m, f_g_m, exflag_g_m, output_g_m] = fmincon(@(ASi)ObjFun(ASi, par, 'multiple_4'), AS0_m, [], [], [], [], ...
                                                    lb_m, ub_m, @(ASi)constraints(ASi, par, 'multiple_4'), options_G_m);
% Display the solution:
fprintf('-------------------------Multiple shooting method with gradient---------------------------------------\n');
disp('The optimised initial state is:')
fprintf('x = %.6f [DU] \n', AS0_opt_g_m(1));
fprintf('y = %.6f [DU] \n', AS0_opt_g_m(2));
fprintf('vx = %.6f [VU] \n', AS0_opt_g_m(3));
fprintf('vy = %.6f [VU] \n', AS0_opt_g_m(4));
fprintf('The initial time is: %.3f [TU] \n', AS0_opt_g_m(end-1));
fprintf('The final time is: %.3f [TU] \n', AS0_opt_g_m(end));
fprintf('The corresponding time of flight is: %.3f \n', AS0_opt_g_m(end)-AS0_opt_g_m(end-1));  % Display the time of flight
fprintf('The optimised value of the objective function (i.e. the cost) is: %.4f \n', f_g_m);  % Display the objective function value
fprintf('Quantities expressed @EMB Earth-Moon rotating frame. \n \n \n ');
 
% Propagation and plot of the found solution:
t_vect_sol = [AS0_opt_g_m(17,1); zeros(3,1)];


% Plot the first node
figure

% Plot the first node
plot(AS0_opt_g_m(1), AS0_opt_g_m(2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, ...
    'DisplayName', 'Nodes')  % Plot the first node as a red dot
hold on

% Plot the trajectories and subsequent nodes
for i = 1:N-1
    t_vect_sol(i+1) = AS0_opt_g_m(17,1) + (i-1+1)/(N-1)*(AS0_opt_g_m(18,1) - AS0_opt_g_m(17,1));  % Time vector for the solution trajectory
    [Sf_g_m, ~, ~, xx_g_m, tt_g_m] = propagation(t_vect_sol(i), AS0_opt_g_m(4*(i-1)+1:4*i,1), t_vect_sol(i+1), par, 'xyPBRFBP_STM_ROT');
    
    if i==1
    % Plot the i-th trajectory
    plot(xx_g_m(:,1), xx_g_m(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410], ...
        'DisplayName', 'Trajectory segments')  % Plot the trajectory for the first segment
    else
            % Plot the i-th trajectory
    plot(xx_g_m(:,1), xx_g_m(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410], ...
        'HandleVisibility', 'off')  % Plot subsequent trajectories without adding to legend
    end
    % Plot the (i+1)-th node
    plot(Sf_g_m(1), Sf_g_m(2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, ...
        'HandleVisibility', 'off')  % Plot the subsequent nodes
end

 
% Add Moon and Earth to the plot
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 12, ...
    'DisplayName', 'Moon')  % Plot the Moon
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 12, ...
    'DisplayName', 'Earth')  % Plot the Earth

% Add grid and axes labels
grid on
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % X-axis label
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % Y-axis label

% Adjust the legend for clarity
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10)  % Legend for the plot

% Ensure the aspect ratio is equal for scientific correctness
axis equal  % Maintain equal scaling on both axes
title('Trajectory in @EMB Earth-Moon rotating frame optimized with gradient and multiple shooting', ...
   'Interpreter', 'latex', 'FontSize', 9) 
hold off

% Define circular orbits
theta = linspace(0, 2*pi, 500);  % Define angle for circular orbits

% Near Earth
% Circular orbit around Earth
x_orbit_earth = -par.mu + par.ri * cos(theta); % Centered at (-par.mu, 0)
y_orbit_earth = par.ri * sin(theta);
figure
plot(xx_g_s(:,1), xx_g_s(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Trajectory');  % Plot the trajectory
hold on  % Hold the plot to overlay the next elements
plot(x_orbit_earth, y_orbit_earth, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2, 'DisplayName','Earth Orbit');  % Plot Earth orbit
plot(-par.mu, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8, 'DisplayName', 'Earth');  % Plot Earth position

% Plot the first node near Earth (coincident with the start of the trajectory)
plot(xx_g_s(1,1), xx_g_s(1,2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Node Near Earth');  % First node

% Add subsequent nodes near Earth
for i = 1:N-1
    plot(Sf_g_m(1), Sf_g_m(2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'HandleVisibility', 'off')  % Plot subsequent nodes near Earth
end
axis equal  % Set equal scaling for x and y axes
xlim([-1.2*par.ri, 1.2*par.ri] - par.mu);  % Limit x-axis range around Earth orbit
ylim([-1.2*par.ri, 1.2*par.ri]);  % Limit y-axis range around Earth orbit
legend('show', 'Interpreter', 'latex', 'Location','best', 'FontSize', 10);
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % X-axis label
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % Y-axis label
grid on;  % Turn on the grid
box on;  % Add box around the plot
title('Trajectory in @EMB Earth-Moon rotating frame optimized with gradient and multiple shooting', ...
   'Interpreter', 'latex', 'FontSize', 9) 
set(gca, 'xtick', [], 'ytick', []);  % Remove x and y ticks

% Near Moon
% Circular orbit around Moon
x_orbit_moon = 1 - par.mu + par.rf * cos(theta); % Centered at (1-par.mu, 0)
y_orbit_moon = par.rf * sin(theta);
figure
plot(xx_g_s(:,1), xx_g_s(:,2), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Trajectory');  % Plot the trajectory
hold on  % Hold the plot to overlay the next elements
plot(x_orbit_moon, y_orbit_moon, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2, 'DisplayName','Moon Orbit');  % Plot Moon orbit
plot(1-par.mu, 0, 'o', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 8, 'DisplayName', 'Moon');  % Plot Moon position

% Plot the first node near Moon (coincident with the start of the trajectory)
plot(xx_g_s(1,1), xx_g_s(1,2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Node Near Moon');  % First node

% Add subsequent nodes near Moon
for i = 1:N-1
    plot(Sf_g_m(1), Sf_g_m(2), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'HandleVisibility', 'off')  % Plot subsequent nodes near Moon
end
axis equal  % Set equal scaling for x and y axes
xlim([-1.2*par.rf, 1.2*par.rf] + (1 - par.mu));  % Limit x-axis range around Moon orbit
ylim([-1.2*par.rf, 1.2*par.rf]);  % Limit y-axis range around Moon orbit
grid on;  % Turn on the grid
box on;  % Add box around the plot
xlabel('$x$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % X-axis label
ylabel('$y$ [DU]', 'Interpreter', 'latex', 'FontSize', 16)  % Y-axis label
legend('show', 'Interpreter', 'latex', 'Location','best', 'FontSize', 10);
set(gca, 'xtick', [], 'ytick', []);  % Remove x and y ticks
title('Trajectory in @EMB Earth-Moon rotating frame optimized with gradient and multiple shooting', ...
   'Interpreter', 'latex', 'FontSize', 9) 
hold off  % End the plot hold


% Characterise the final orbit:
[orbitType_g_m, captureType_g_m] = characterise_orbits(Sf_g_m(:,:,end), par);
disp('The achieved orbit can be characterised as follows:') % Display the orbit and capture types
fprintf(' - %s \n', orbitType_g_m); 
fprintf(' - %s \n \n', captureType_g_m);  



%% Exercise 2.4
% Set useful parameters:
frame = 'ECLIPJ2000';

% Initial Moon-EMB-Sun angle:
theta_i = wrapTo2Pi(par.ws*AS0_opt_g_s(5));

% Starting epoch:
t_et_0_TDB = '2024 Sep 28 00:00:00.000 TDB';
t_et_0 = cspice_str2et(t_et_0_TDB);              % Convert starting epoch in ephemeris tim

% Plot the evolution of theta from t_et_0_TDB to epoch_f, to properly set
% the boundaries for the fzero:
epoch_0 = '2024 Sep 28 00:00:00.000';  % Initial epoch
epoch_f = '2024 Dec 28 00:00:00.000';  % Final epoch
t_f = cspice_str2et(epoch_f);          % Set final month for display of the solution
t_vect = linspace(t_et_0, t_f, 1000);  % Time vector
theta = zeros(length(t_vect), 1);      % Pre-allocate theta vector
t_vect_date = cell(length(t_vect), 1); % Pre-allocate cell array for date strings

% Calculate theta and corresponding dates:
for i = 1:length(t_vect)
    [~, theta(i)] = thetafind(t_vect(i), theta_i);
    t_vect_date{i} = cspice_timout(t_vect(i), 'YYYY MON DD HR:MN:SC.####::TDB');
end 

% Convert date strings to datetime format for plotting
date_epochs = datetime(t_vect_date, 'InputFormat','yyyy MMM dd HH:mm:ss.SSSS', 'Format', 'd MMM uuuu');

% Create the plot:
figure;

% Plot the theta evolution with an elegant line style and color
theta_plot = plot(date_epochs, theta, 'b-', 'LineWidth', 1.6, 'Color', [0 0.4470 0.7410], 'DisplayName', '$\theta$ Evolution');
hold on;

% Reference line for theta_i
theta_i_line = yline(theta_i, '--r', 'LineWidth', 1.2, 'Interpreter', 'latex', ...
    'Label', '$\theta_i$', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'Color', [1 0 0]);

% Define vertical lines for 30-day intervals
start_date = datetime(epoch_0, 'InputFormat', 'yyyy MMM dd HH:mm:ss.SSSS');
end_date = datetime(epoch_f, 'InputFormat', 'yyyy MMM dd HH:mm:ss.SSSS');
day_interval = days(30); 
xlines = gobjects(1, length(start_date:day_interval:end_date)); % Pre-allocate for legend
idx = 1;

for d = start_date:day_interval:end_date
    if idx == 1
        % Add the first 30-day interval line to the legend
        xlines(idx) = xline(d, '--k', 'LineWidth', 0.8, 'Alpha', 0.4, 'DisplayName', '30-day Intervals');
    else
        % Other lines are plotted without adding to the legend
        xlines(idx) = xline(d, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);
    end
    idx = idx + 1;
end

% Add title and labels
title('Sun-EMB Angle Evolution', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
xlabel('Date [TDB]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\theta$ [rad]', 'FontSize', 14, 'Interpreter', 'latex');

% Customize x-axis ticks
xticks = linspace(date_epochs(1), date_epochs(end), 6); % Divide x-axis into six intervals
xticklabels_cell = arrayfun(@(d) sprintf('$%s$', datetime(d, 'InputFormat', 'yyyy MMM dd HH:mm:ss.SSSS')), ...
    xticks, 'UniformOutput', false);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels_cell, 'TickLabelInterpreter', 'latex');

% Add legend with improved formatting
legend([theta_plot, theta_i_line, xlines(1)], {'$\theta$ Evolution', '$\theta_i$ Reference', '30-day intervals'}, ...
       'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');

% Refine general plot aesthetics
set(gca, 'FontSize', 12, 'GridLineStyle', '-', 'LineWidth', 1);
hold off;
% From the plot, it can be note that theta values cross theta_i value in a
% date in between 8 Oct and 16 Oct. Those bounds are passed to fzero
% algorithm to find a solution, since are the one in which the function
% changes its sign.

% Solve for initial time matching the specified angle and in the retrieved 
% epoches:
options_fzero = optimset('TolX', 1e-15);
t_et_i = fzero(@(t) thetafind(t, theta_i), [cspice_str2et('2024 Oct 8 00:00:00.000 TDB'), cspice_str2et('2024 Oct 16 00:00:00.000 TDB')], options_fzero);
t_et_i_TDB = cspice_timout(t_et_i, 'YYYY-MM-DD HR:MN:SC.###::TDB');
t_et_i_UTC = cspice_timout(t_et_i, 'YYYY-MM-DD HR:MN:SC.###::UTC');


% Define final epoch time based on initial and scaled final time:
tof = (AS0_opt_g_s(6) - AS0_opt_g_s(5))*par.TU*24*3600;
t_et_f = t_et_i + tof;                                                      % Final time in ET
t_et_f_TDB = cspice_timout(t_et_f, 'YYYY-MM-DD HR:MN:SC.###::TDB');         % Convert final time to TDB format
t_et_f_UTC = cspice_timout(t_et_f, 'YYYY-MM-DD HR:MN:SC.###::UTC');         % Convert final time to TDB format


% Print the initial and final epoch times in TDB format
fprintf('\n------------------------------Epoch finding-------------------------------- \n')
fprintf('Initial epoch: %s (UTC)\n ', t_et_i_UTC);
fprintf('Final epoch: %s(UTC)\n \n', t_et_f_UTC);


% N-body propagation:
fprintf('------------------------------N-body propagation--------------------------------------\n')

% Propagate the orbit firstly directly in ECLIPJ2000 coordinates and in rotational frame
% Define celestial bodies for n-body problem:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data:
bodies = nbody_init(labels);

% Initial state vector transformation (rotating frame to ECI):
AS0_4_ECI = zeros(6, 1);
AS0_4_ECI([1,2,4,5]) = ROT2ECI(AS0_opt_g_s(1:4)', AS0_opt_g_s(5), par); % Apply reference frame transformation
AS0_4_ECI(1:3) = par.DU * AS0_4_ECI(1:3); % Scale position
AS0_4_ECI(4:6) = par.VU * AS0_4_ECI(4:6); % Scale velocity

% Print position vector
fprintf('Position Vector (ECI): [%.8f, %.8f, %.8f] km\n', AS0_4_ECI(1), AS0_4_ECI(2), AS0_4_ECI(3));

% Print velocity vector
fprintf('Velocity Vector (ECI): [%.8f, %.8f, %.8f] km/s\n \n', AS0_4_ECI(4), AS0_4_ECI(5), AS0_4_ECI(6));

% Propagation in ECLIPJ2000 coordinates:
[~, ~, ~, xx_ECI, tt_plot] = propagation(t_et_i, AS0_4_ECI, t_et_f, par, 'xyz_n_body', bodies, frame, 'Earth');

% Convert rotated trajectory to ECI:
xx_rot = ROT2ECI(xx_g_s, tt_g_s, par);

% Moon trajectory in PBRFBP calculation:
t_moon = linspace(AS0_opt_g_s(5), AS0_opt_g_s(6), size(xx_rot, 1)); % Moon time vector
xx_moon_ROT = [(1-par.mu) * ones(size(xx_rot, 1), 1), zeros(size(xx_rot, 1), 3)];
xx_moon = ROT2ECI(xx_moon_ROT, t_moon, par);

% Convert propagation time to absolute Ephemeris Time:
t_moon_kernel = (tt_g_s - tt_g_s(1))*par.TU*24*3600 + t_et_i;

% Moon trajectory extraction from the kernels (in ECI):
xx_moon_kernel = cspice_spkezr('Moon', t_moon_kernel', 'ECLIPJ2000', 'none', 'Earth');

% Plotting the results in ECI:
figure;
hold on;
theta = linspace(0, 2*pi, 250);
plot3(par.DU * cos(theta), par.DU* sin(theta), zeros(size(theta)), ...
    '--', 'LineWidth', 1, 'DisplayName', 'Moon''s Trajectory PBRFBP', 'Color', [1 0 0]);
plot3(xx_moon_kernel(1, :), xx_moon_kernel(2, :), xx_moon_kernel(3, :), ...
      'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Moon real trajectory (in ECLIPJ2000)');
plot3(xx_ECI(:, 1), xx_ECI(:, 2), xx_ECI(:, 3), ...
      'Color', [0.6 0 0.6], 'LineWidth', 1.5, 'DisplayName', 'N-Body propagated trajectory');
plot3(xx_rot(:,1) * par.DU, xx_rot(:,2) * par.DU, zeros(size(xx_rot, 1), 1), ...
      'Color', [0.00,0.45,0.74], 'LineWidth', 1.5, 'DisplayName', '2D PBRFBP propagated trajectory');
plot3(0, 0, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8, 'DisplayName', 'Earth')
axis equal;
grid on;
xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ [km]', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$z$ [km]', 'Interpreter', 'latex', 'FontSize', 14);
title('Spacecraft Trajectories in ECI', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best', 'Interpreter', 'latex');
view(3);
hold off;

% Dimensionalise:
AAS = [(xx_g_s(:, 1) + par.mu) * par.DU, xx_g_s(:, 2) * par.DU, zeros(length(xx_g_s), 1), ...
     xx_g_s(:, 3) * par.VU, xx_g_s(:, 4) * par.VU, zeros(length(xx_g_s), 1)];

% Moon trajectory extraction from the kernels (in ECLIPJ2000):
xx_moon_kernel = cspice_spkezr('Moon', t_moon_kernel', 'ECLIPJ2000', 'none', 'Earth');

% Initialize transformed trajectory
xx_obj_ECLIPJ2000 = zeros(6, length(t_moon_kernel));

% Transform trajectory from rotating frame to J2000
for i = 1:length(t_moon_kernel)
    r_moon_kernel = xx_moon_kernel(1:3, i);
    v_moon_kernel = xx_moon_kernel(4:6, i);

    % Compute rotating frame axes
    X_tilde = r_moon_kernel / norm(r_moon_kernel);
    Z_tilde = cross(r_moon_kernel, v_moon_kernel) / norm(cross(r_moon_kernel, v_moon_kernel));
    Y_tilde = cross(Z_tilde, X_tilde);
    R_RT2ECLIPJ2000 = [X_tilde, Y_tilde, Z_tilde];

    % Transform position and velocity
    theta_dot = norm(cross(r_moon_kernel, v_moon_kernel)) / norm(r_moon_kernel)^2;
    Rv_ECLIPJ2000 = [theta_dot * R_RT2ECLIPJ2000(1, 2), -theta_dot * R_RT2ECLIPJ2000(1, 1), 0, R_RT2ECLIPJ2000(1, :); ...
                theta_dot * R_RT2ECLIPJ2000(2, 2), -theta_dot * R_RT2ECLIPJ2000(2, 1), 0, R_RT2ECLIPJ2000(2, :); ...
                theta_dot * R_RT2ECLIPJ2000(3, 2), -theta_dot * R_RT2ECLIPJ2000(3, 1), 0, R_RT2ECLIPJ2000(3, :)];
    xx_obj_ECLIPJ2000(:, i) = [R_RT2ECLIPJ2000 zeros(3); Rv_ECLIPJ2000] * AAS(i, :)';
end

% N-body propagation
[~, ~, ~, xx_nbody_tr, tt_nbody_tr] = propagation(t_et_i, xx_obj_ECLIPJ2000(:, 1), t_et_f, par, 'xyz_n_body', bodies, frame, 'Earth');

% Plot:
figure; hold on;
plot3(0, 0, 0, 'o', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 8, 'DisplayName', 'Earth')
plot3(par.DU * cos(theta), par.DU * sin(theta), zeros(size(theta)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Moon Trajectory (2D PBRFBP)');
plot3(xx_nbody_tr(:, 1), xx_nbody_tr(:, 2), xx_nbody_tr(:, 3), 'Color', [0.6 0 0.6], 'LineWidth', 1.5, 'DisplayName', 'N-Body propagated trajectory');
plot3(xx_obj_ECLIPJ2000(1, :), xx_obj_ECLIPJ2000(2, :), xx_obj_ECLIPJ2000(3, :), 'Color', [0.00,0.45,0.74], 'LineWidth', 1.5, 'DisplayName', '2D PBRFBP propagated trajectory');
plot3(xx_moon_kernel(1, :), xx_moon_kernel(2, :), xx_moon_kernel(3, :), 'DisplayName', 'Moon real trajectory (in ECLIPJ2000)', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('x [km]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('y [km]', 'FontSize', 12, 'Interpreter', 'latex');
zlabel('z [km]', 'FontSize', 12, 'Interpreter', 'latex');
title('Spacecraft Trajectories @Earth ECLIPTICJ2000', 'FontSize', 14);
legend('show', 'Location', 'southeast' , 'FontSize', 10, 'Interpreter', 'latex');
grid on;
axis equal
view(3);


%% Functions

% -------------------------------------------------------------------------
function [AS_first_guess] = first_guess(alpha, beta, delta, ti0, par)
% DESCRIPTION:
% Finds the first guess solution for the Sun-(Earth+Moon) Planar Bicircular 
% Restricted Four Body Problem for a fixed ri problem, using four input parameters.
% 
% PROTOTYPE:
% [AS_first_guess] = first_guess(alpha, beta, delta, ti0, par)
%
% INPUT:
% alpha [1x1]     Angle on the Earth circular parking orbit of radius ri  [rad]
% beta [1x1]      Initial-to-circular velocity ratio                      [-]
% ti0  [1x1]      Initial time                                            [s]
% delta [1x1]     Transfer duration                                       [s]
% par   [struct]  Structure containing constants (mu, ri, etc.)
%
% OUTPUT:
% AS_first_guess [8x1]  Augmented state of the first guess solution
%                       (initial position, velocity, and times)

    % Rename variables:
    mu = par.mu;         % Gravitational parameter
    r0 = par.ri;         % Radius of Earth parking orbit
    v0 = beta*sqrt((1-mu)/r0); % Initial velocity

    % Calculation of the initial guess for initial conditions:
    x0 = r0*cos(alpha) - mu;   % Initial x-coordinate
    y0 = r0*sin(alpha);        % Initial y-coordinate
    dx0 = -(v0 - r0)*sin(alpha); % Initial velocity in x-direction
    dy0 = (v0 - r0)*cos(alpha);  % Initial velocity in y-direction
    S0 = [x0; y0; dx0; dy0];    % State vector

    % The final time is:
    tf0 = ti0 + delta;          % Final time

    % The augmented state is:
    AS_first_guess = [S0; ti0; tf0]; % Augmented state output

end

% -------------------------------------------------------------------------
function [dSdt] = xyPBRFBP_STM_ROT(t, S0, par)
% DESCRIPTION:
% Computes the right-hand side of the state vector dynamics for the Planar 
% Bicircular Restricted Four Body Problem, considering the State Transition 
% Matrix (STM) in a rotating frame.
%
% PROTOTYPE:
% [dSdt] = xyPBRFBP_STM_ROT(t, S0, par)
%
% INPUT:
% S0 [20x1]  Initial state vector (position, velocity, and STM)
% par [struct]  Structure containing parameters like mu, ms, ws, rho
%
% OUTPUT:
% dSdt [20x1]  Derivative of the state vector (dynamics)

    
    % Extract variables
    x  = S0(1);   % x-coordinate of position vector
    y  = S0(2);   % y-coordinate of position vector
    vx = S0(3);   % x-coordinate of velocity vector
    vy = S0(4);   % x-coordinate of position vector
    
    % Define preliminarly some useful parameters:
    mu = par.mu;        % Gravitational parameter
    ms = par.ms;        % [-]: Scaled mass of the Sun
    ws = par.ws;        % [-]: Scaled angular velocity of the Sun
    rho = par.rho;      % [-]: Scaled Sun-(Earth+Moon) distance

    % Put PHI in matrix form
    Phi = reshape(S0(5:end),4,4);

    % Calculate the important radii (same as before, but now reusable):
    r1 = (mu + x)^2 + y^2;
    r2 = (mu + x - 1)^2 + y^2;
    r3 = (x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2;
    
    % Compute r1^(3/2), r2^(3/2), r3^(3/2) for more consistent terms
    r1_3_2 = r1^(3/2);  % r1^(3/2)
    r2_3_2 = r2^(3/2);  % r2^(3/2)
    r3_3_2 = r3^(3/2);  % r3^(3/2)
    
    % Compute the derivative of the potential consistently
    dOM4dx = x - (mu*(mu + x - 1))/(r2_3_2) - (ms*cos(t*ws))/rho^2 + ((mu + x)*(mu - 1))/(r1_3_2) - (ms*(x - rho*cos(t*ws)))/r3_3_2;
    dOM4dy = y - (ms*sin(t*ws))/rho^2 - (mu*y)/r2_3_2 - (ms*(y - rho*sin(t*ws)))/r3_3_2 + (y*(mu - 1))/r1_3_2;
    
    % Assemble the matrix A(t)=dfdx 4x4 matrix
    dfdx = [                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                             0,  1, 0;
                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                             0,  0, 1;
           (mu - 1)/((mu + x)^2 + y^2)^(3/2) - ms/((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(t*ws))^2)/((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(5/2) + 1,                                                                      (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) + (3*ms*(2*x - 2*rho*cos(t*ws))*(2*y - 2*rho*sin(t*ws)))/(4*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)),  0, 2;
                                                                                                  (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) + (3*ms*(2*x - 2*rho*cos(t*ws))*(2*y - 2*rho*sin(t*ws)))/(4*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - ms/((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(t*ws))^2)/((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, -2, 0];
 
    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

    % Assemble right-hand side
    dSdt = zeros(20,1);

    % Position and velocity derivatives
    dSdt(1:2) = S0(3:4);          % Position derivatives (velocity)
    dSdt(3)   = dOM4dx + 2*vy;    % Velocity x-derivative with Coriolis term
    dSdt(4)   = dOM4dy - 2*vx;    % Velocity y-derivative with Coriolis term

    % State Transition Matrix (STM) derivatives
    dSdt(5:end) = Phidot(:);  % Flatten and assign STM derivatives

    
end

% -------------------------------------------------------------------------
function [Sf, PHIf, tf, xx, tt]  = propagation(t0, S0, tf, par, label, bodies, frame, center)
% [xf, PHIf, tf, xx, tt]  = propagation(t0, S0, tf, par, label)
% Description:
% Perform the propagation in a reference system specified and allowed by
% "label". This code uses ode78 for the propagation, since it best suites
% the astrodynamical problem.
% 
% Input:
% t0: [1x1] initial time for the propagation [s]
% S0: [??x1] initial state for the propagation. If label='xyPBRFBP_STM_ROT'
%     is chosen, size(S0)=[4x1]. If label='xyz_n_body' is chosen, the 
%     initial state is necessary in 3D, thus size(S0)=[6x1]
% tf: [1x1] Final time for the propagation
% par: [struct] Structure containing all the useful quantities to perform the
%            computation
% label: [str] Choose between 'xyPBRFBP_STM_ROT' and 'xyz_n_body'. If
%              'xyz_n_body' is chosen, then specify the other following 
%              inputs.
% bodies : [1,n] cell-array with struct elements containing the following
%                fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
% 
% frame: [str] reference frame in which the propagation must be performed
% centre: [str] centre with respect to which the propagation must be
%         performed
%
% Outputs:
% Sf: [6x1] final state at tf
% PHIf: [6x6] STM at tf. This output is set to null matrix in case 
%       label='xyz_n_body' is chosen
% tf: [1x1] final time 
% xx: [6x?)] matrix containing the evolution of the state vector along the
%     integration
% tt: [1x??] vector containing the evolution of the time index along the
%     propagaion

    if strcmp(label, 'xyPBRFBP_STM_ROT') || nargin < 6
    
        % Initialize State Transition Matrix at t0
        Phi0 = eye(4);
    
        % Append to initial conditions the conditions for the STM
        x0Phi0 = [S0; Phi0(:)];
        
        % Perform integration
        options = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,2), 1e-12*ones(1,2), 1e-12*ones(1,16)]);
        % options = odeset('reltol', 1e-12, 'abstol', 1e-12);

        [tt, xx] = ode78(@(t,x) xyPBRFBP_STM_ROT(t, x, par), [t0 tf], x0Phi0, options);
    
        % Extract state vector and State Transition Matrix
        Sf = xx(end,1:4)';
        PHIf = reshape(xx(end,5:end),4,4);
        tf = tt(end);

    elseif strcmp(label, 'xyz_n_body') || nargin >= 6
        
        % Perform integration
        options = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,3), 1e-12*ones(1,3)]);
        % options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [tt, xx] = ode78(@(t,x) nbody_shift_rhs(t, x, bodies, frame, center), [t0 tf], S0, options);
    
        % Extract (final) state vector and State Transition Matrix:
        Sf = xx(end,1:6)';
        PHIf = zeros(6,6);
        tf = tt(end);

    else 
        error('The label chosen is not valid')
    end 

end

% -------------------------------------------------------------------------
function S_EC = ROT2ECI(S_ROT, t, par, str)
% DESCRIPTION:
% Converts the state vector from the rotating frame to the Earth-centered 
% inertial (ECI) reference frame.
%
% PROTOTYPE:
% S_EC = ROT2ECI(S_ROT, t, par, str)
%
% INPUT:
% S_ROT [Nx4]      State vector in the rotating frame                [x, y, vx, vy] (km, km/s)
% t    [Nx1]       Time vector for transformation                    [s]
% par  [struct]    Structure containing physical parameters:
%                  mu [1x1]  Gravitational parameter                [km^3/s^2]
% str  [string]    Optional string ('moon' for Moon-specific transformation)
%
% OUTPUT:
% S_EC [Nx4]       State vector in the Earth-Centered Inertial (ECI) frame
%                  [X, Y, VX, VY] (km, km/s)
%

    % Extract the gravitational parameter from the input structure
    mu = par.mu;
    
    % Initialize the output state vector
    S_EC = zeros(size(S_ROT));

    % Loop through each time step to perform the conversion
    for i=1:length(t)
        
        % Allocate the variables for position and velocity components
        x = S_ROT(i, 1);   % Position in x-coordinate in rotating frame
        y = S_ROT(i, 2);   % Position in y-coordinate in rotating frame
        vx = S_ROT(i, 3);  % Velocity in x-coordinate in rotating frame
        vy = S_ROT(i, 4);  % Velocity in y-coordinate in rotating frame
        
        % Perform the transformation from rotating frame to ECI frame
        X = (x + mu)*cos(t(i)) - y*sin(t(i));  % Transformed x-coordinate
        Y = (x + mu)*sin(t(i)) + y*cos(t(i));  % Transformed y-coordinate
        VX = (vx - y)*cos(t(i)) - (vy + x + mu)*sin(t(i));  % Transformed velocity in x
        VY = (vx - y)*sin(t(i)) + (vy + x + mu)*cos(t(i));  % Transformed velocity in y
        
        % If the 'moon' option is specified, modify the transformation
        if nargin > 3 && strcmp(str, 'moon')
            X = (x + mu)*cos(t(i)) + y*sin(t(i));  % Adjusted x for Moon
            Y = (x + mu)*(-sin(t(i))) + y*cos(t(i));  % Adjusted y for Moon
        end
        
        % Store the results in the output state vector
        S_EC(i, 1:4) = [X; Y; VX; VY];
    end
end

% -------------------------------------------------------------------------
function [F,G] = ObjFun(ASi, par, shooting_technique)
% DESCRIPTION:
% This function defines the objective function and its gradient (if requested) 
% for the Sun-(Earth+Moon) Planar Bicircular Restricted Four-Body Problem (PBRFBP).
% It calculates the Δv at both the initial and final states based on the shooting 
% technique (simple or multiple_4) and computes the gradient if requested.
% 
% PROTOTYPE:
% [F, G] = ObjFun(ASi, par, shooting_technique)
%
% INPUTS:
% ASi [1x6] or [1x18] array: Variables of the optimization problem. 
%      If shooting_technique == 'simple', ASi must be [6x1] or [1x6].
%      If shooting_technique == 'multiple_4', ASi must be [18x1] or [1x18].
% par [struc]: Structure containing the parameters for the model.
% shooting_technique [str]: The method to use, either 'simple' or 'multiple_4'.
%
% OUTPUTS:
% F [1x1]: Objective function value (sum of initial and final Δv).
% G [1x6] or [1x18]: Gradient of the objective function with respect to the state.
%                     Empty if gradient is not requested.

% Define preliminarly some useful parameters:
mu = par.mu;        % Gravitational parameter
ms = par.ms;        % [-]: Scaled mass of the Sun
ws = par.ws;        % [-]: Scaled angular velocity of the Sun
rho = par.rho;      % [-]: Scaled Sun-(Earth+Moon) distance
ri = par.ri;        % [-]: Initial parking orbit radius
rf = par.rf;        % [-]: Final parking orbit radius
   
    if strcmp(shooting_technique, 'simple') == 1 % Simple shooting technique
    
    % Allocate the variables conveniently:
    xi = ASi(1);                   % x position coordinate of initial state
    yi = ASi(2);                   % y position coordinate of initial state
    vxi = ASi(3);                  % x velocity coordinate of initial state
    vyi = ASi(4);                  % y velocity coordinate of initial state
    S0 = [xi, yi, vxi, vyi]';      % Initial state
    ti = ASi(5);                   % Initial time
    tf = ASi(6);                   % Final time

    % Calculate the final state and allocate the relative variables:
    [Sf, PHIf, tf, ~, ~] = propagation(ti, S0, tf, par, 'xyPBRFBP_STM_ROT');
    xf = Sf(1);                    % x position coordinate of final state
    yf = Sf(2);                    % y position coordinate of final state
    vxf = Sf(3);                   % x velocity coordinate of final state
    vyf = Sf(4);                   % y velocity coordinate of final state

    % Calculate the initial and final impulses:
    Dvi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - sqrt((1 - mu) / ri);
    Dvf = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);

    % The objective function:
    F = Dvi + Dvf;

    if nargout > 1 % gradient calculations requested
        
        % Allocate space for the gradient
        G = zeros(1, 6);

        % Define the common expressions for gradient calculation:
        r1_sq = (vxi - yi)^2 + (mu + vyi + xi)^2;
        r1_inv = 1 / sqrt(r1_sq);  % Reciprocal of the square root for Dvi

        r2_sq = (vxf - yf)^2 + (mu + vyf + xf - 1)^2;
        r2_inv = 1 / sqrt(r2_sq);  % Reciprocal of the square root for Dvf

        % Define the derivative of the impulses with respect to the initial
        % adn final state rspectively:
        dDvi_dSi = [(mu + vyi + xi), -(vxi - yi), (vxi - yi), (mu + vyi + xi)]*r1_inv;
        dDvf_dSf = [(mu + vyf + xf - 1), -(vxf - yf), (vxf - yf), (mu + vyf + xf - 1)]*r2_inv;

        % Compute the dynamics vector at initial and final state:
        dOM4dx_i = xi - (mu * (2 * mu + 2 * xi - 2)) / (2 * ((mu + xi - 1)^2 + yi^2)^(3 / 2)) - (ms * cos(ti * ws)) / rho^2 + ((2 * mu + 2 * xi) * (mu - 1)) / (2 * ((mu + xi)^2 + yi^2)^(3 / 2)) - (ms * (2 * xi - 2 * rho * cos(ti * ws))) / (2 * ((xi - rho * cos(ti * ws))^2 + (yi - rho * sin(ti * ws))^2)^(3 / 2));
        dOM4dy_i = yi - (ms * sin(ti * ws)) / rho^2 - (mu * yi) / ((mu + xi - 1)^2 + yi^2)^(3 / 2) - (ms * (2 * yi - 2 * rho * sin(ti * ws))) / (2 * ((xi - rho * cos(ti * ws))^2 + (yi - rho * sin(ti * ws))^2)^(3 / 2)) + (yi * (mu - 1)) / ((mu + xi)^2 + yi^2)^(3 / 2);
        dOM4dx_f = xf - (mu * (2 * mu + 2 * xf - 2)) / (2 * ((mu + xf - 1)^2 + yf^2)^(3 / 2)) - (ms * cos(tf * ws)) / rho^2 + ((2 * mu + 2 * xf) * (mu - 1)) / (2 * ((mu + xf)^2 + yf^2)^(3 / 2)) - (ms * (2 * xf - 2 * rho * cos(tf * ws))) / (2 * ((xf - rho * cos(tf * ws))^2 + (yf - rho * sin(tf * ws))^2)^(3 / 2));
        dOM4dy_f = yf - (ms * sin(tf * ws)) / rho^2 - (mu * yf) / ((mu + xf - 1)^2 + yf^2)^(3 / 2) - (ms * (2 * yf - 2 * rho * sin(tf * ws))) / (2 * ((xf - rho * cos(tf * ws))^2 + (yf - rho * sin(tf * ws))^2)^(3 / 2)) + (yf * (mu - 1)) / ((mu + xf)^2 + yf^2)^(3 / 2);
        f_i = [vxi; vyi; dOM4dx_i + 2*vyi; dOM4dy_i - 2*vxi];
        f_f = [vxf; vyf; dOM4dx_f + 2*vyf; dOM4dy_f - 2*vxf];
        
        % Calculate the gradient of the function:
        G(1:4) = dDvi_dSi + dDvf_dSf*PHIf;
        G(5) = -dDvf_dSf*PHIf*f_i;
        G(6) = dDvf_dSf*f_f;

     end
        
    
    elseif strcmp(shooting_technique, 'multiple_4') % Multiple shooting techinque with N=4
        
        % Allocate conveniently the variables:
        Si = ASi(1:4);   % Initial state
        xi = Si(1);      % x position coordinate of initial state
        yi = Si(2);      % y position coordinate of initial state
        vxi = Si(3);     % x velocity coordinate of initial state
        vyi = Si(4);     % y position coordinate of initial state
        Sf = ASi(13:16); % Final state
        xf = Sf(1);      % x position coordinate of final state
        yf = Sf(2);      % y position coordinate of final state
        vxf = Sf(3);     % x velocity coordinate of final state
        vyf = Sf(4);     % y position coordinate of final state

        
        % Calculate the initial and final impulses:
        Dvi = sqrt((vxi-yi)^2+(vyi+xi+mu)^2) - sqrt((1-mu)/ri);
        Dvf = sqrt((vxf-yf)^2+(vyf+xf+mu-1)^2) - sqrt(mu/rf);
    
        % The objective function is:
        F = Dvi + Dvf;

        if nargout>1 % gradient calculations requested

            % Allocate the variable:
            G = zeros(1,18);

            % From the reference it is possible to compute that:
            G(1:4) = [vyi+xi+mu, yi-vxi, vxi-yi, vyi+xi+mu]./(sqrt((vxi-yi)^2+(vyi+xi+mu)^2));
            G(5:12) = zeros(1,8);
            G(13:16) = [vyf+xf+mu-1, yf-vxf, vxf-yf, vyf+xf+mu-1]./(sqrt((vxf-yf)^2+(vyf+xf+mu-1)^2));
            G(17:18) = zeros(1,2);
 
        end 

    else
        error('Unknown method specified. Choose between ''simple'' and ''multiple_4''')
    end 


end

% -------------------------------------------------------------------------
function [c, ceq, gc, gceq] = constraints(ASi, par, shooting_technique)
% DESCRIPTION:
% This function defines the objective function and its gradient (if requested) 
% for the Sun-(Earth+Moon) Planar Bicircular Restricted Four-Body Problem (PBRFBP).
% It calculates the Δv at both the initial and final states based on the shooting 
% technique (simple or multiple_4) and computes the gradient if requested.
% 
% PROTOTYPE:
% [F, G] = ObjFun(ASi, par, shooting_technique)
%
% INPUTS:
% ASi [1x6] or [1x18] array: Variables of the optimization problem. 
%      If shooting_technique == 'simple', ASi must be [6x1] or [1x6].
%      If shooting_technique == 'multiple_4', ASi must be [18x1] or [1x18].
% par [struc]: Structure containing the parameters for the model.
% shooting_technique [str]: The method to use, either 'simple' or 'multiple_4'.
% OUTPUTS:
% F [1x1]: Objective function value (sum of initial and final Δv).
% G [1x6] or [1x18]: Gradient of the objective function with respect to the state.
%                     Empty if gradient is not requested.

% Define preliminarly some useful parameters:
mu = par.mu;           % Gravitational parameter
ms = par.ms;           % [-]: Scaled mass of the Sun
ws = par.ws;           % [-]: Scaled angular velocity of the Sun
rho = par.rho;         % [-]: Scaled Sun-(Earth+Moon) distance
ri = par.ri;           % [-]: Initial parking orbit radius
rf = par.rf;           % [-]: Final parking orbit radius
E_Radii = par.E_Radii; % [km]: Vector of radii of Earth
M_Radii = par.M_Radii; % [km]: Vector of radii of Moon
DU = par.DU;           % [km]: Distance Unit
   
    if strcmp(shooting_technique, 'simple')==1 % Simple shooting technique

        % Allocate conveniently the variables:
        xi = ASi(1);
        yi = ASi(2);
        vxi = ASi(3);
        vyi = ASi(4);
        S0 = [xi, yi, vxi, vyi]';
        ti = ASi(5);
        tf = ASi(6);
       
        % Calculate the final state:
        [Sf, PHIf, ~, ~, ~]  = propagation(ti, S0, tf, par, 'xyPBRFBP_STM_ROT');
    
        % Allocate conveniently the variables:
        xf = Sf(1);
        yf = Sf(2);
        vxf = Sf(3);
        vyf = Sf(4);
    
        % Inequality constraints (not present):
        c = [];
    
        % Define the initial and final conditions for the described transfer:
        Psi1 = [(xi+mu)^2 + yi^2-ri^2; (xi+mu)*(vxi-yi)+yi*(vyi+xi+mu)];
        Psi2 = [(xf+mu-1)^2+yf^2-rf^2; (xf+mu-1)*(vxf-yf)+yf*(vyf+xf+mu-1)];
               
        % The equality constraints:
        ceq = [Psi1; Psi2];
    
        if nargout>2
    
            % Define and calculate some useful quantities:
            dOM4dx_i = xi - (mu * (2 * mu + 2 * xi - 2)) / (2 * ((mu + xi - 1)^2 + yi^2)^(3 / 2)) - (ms * cos(ti * ws)) / rho^2 + ((2 * mu + 2 * xi) * (mu - 1)) / (2 * ((mu + xi)^2 + yi^2)^(3 / 2)) - (ms * (2 * xi - 2 * rho * cos(ti * ws))) / (2 * ((xi - rho * cos(ti * ws))^2 + (yi - rho * sin(ti * ws))^2)^(3 / 2));
            dOM4dy_i = yi - (ms * sin(ti * ws)) / rho^2 - (mu * yi) / ((mu + xi - 1)^2 + yi^2)^(3 / 2) - (ms * (2 * yi - 2 * rho * sin(ti * ws))) / (2 * ((xi - rho * cos(ti * ws))^2 + (yi - rho * sin(ti * ws))^2)^(3 / 2)) + (yi * (mu - 1)) / ((mu + xi)^2 + yi^2)^(3 / 2);
    
            dOM4dx_f = xf - (mu * (2 * mu + 2 * xf - 2)) / (2 * ((mu + xf - 1)^2 + yf^2)^(3 / 2)) - (ms * cos(tf * ws)) / rho^2 + ((2 * mu + 2 * xf) * (mu - 1)) / (2 * ((mu + xf)^2 + yf^2)^(3 / 2)) - (ms * (2 * xf - 2 * rho * cos(tf * ws))) / (2 * ((xf - rho * cos(tf * ws))^2 + (yf - rho * sin(tf * ws))^2)^(3 / 2));
            dOM4dy_f = yf - (ms * sin(tf * ws)) / rho^2 - (mu * yf) / ((mu + xf - 1)^2 + yf^2)^(3 / 2) - (ms * (2 * yf - 2 * rho * sin(tf * ws))) / (2 * ((xf - rho * cos(tf * ws))^2 + (yf - rho * sin(tf * ws))^2)^(3 / 2)) + (yf * (mu - 1)) / ((mu + xf)^2 + yf^2)^(3 / 2);
    
            % The dynamics at initial and final state:
            f_i = [vxi; vyi; dOM4dx_i + 2*vyi; dOM4dy_i - 2*vxi];
            f_f = [vxf; vyf; dOM4dx_f + 2*vyf; dOM4dy_f - 2*vxf];

            % Inequality constraint gradient (not present):
            gc = [];
    
            % Allocate the space for the gradient of the equality constraints:
            gceq = zeros(4,6);
      
            % Define some useful quantities:
            dceqi_dSi = [2*mu + 2*xi, 2*yi, 0,  0; vxi,  vyi, mu + xi, yi];
            dceqf_dSf = [2*mu + 2*xf - 2, 2*yf, 0,  0; vxf,  vyf, mu + xf - 1, yf];
            
            % Assemble the gradient of the equality constraints:
            gceq(1:4, 1:4) = [dceqi_dSi; dceqf_dSf*PHIf];
            gceq(1:4, 5) = [zeros(2,1); -dceqf_dSf*PHIf*f_i];
            gceq(1:4, 6) = [zeros(2,1); dceqf_dSf*f_f];
            gceq= gceq';
    
        end 
    
    elseif strcmp(shooting_technique, 'multiple_4')

        N = 4; % Numbers of shooting

        % Allocate conveniently the variables:
        Si = ASi(1:4);   % Initial state
        xi = Si(1);      % x position coordinate of initial state
        yi = Si(2);      % y position coordinate of initial state
        vxi = Si(3);     % x velocity coordinate of initial state
        vyi = Si(4);     % y position coordinate of initial state
        ti = ASi(17);    % Initial time
        
        Sf = ASi(13:16); % Final state
        xf = Sf(1);      % x position coordinate of final state
        yf = Sf(2);      % y position coordinate of final state
        vxf = Sf(3);     % x velocity coordinate of final state
        vyf = Sf(4);     % y position coordinate of final state
        tf = ASi(18);    % Final time

         % Find intermediate conditions:
        S2 = ASi(5:8);                 % State at time t2          
        x2 = S2(1);                    % x-coordinate of position at time t2
        y2 = S2(2);                    % y-coordinate of position at time t2
        vx2 = S2(3);                   % x-coordinate of velocity at time t2
        vy2 = S2(4);                   % y-coordinate of velocity at time t2
        t2 = ti + (2-1)/(N-1)*(tf-ti); % Time t2
        S3 = ASi(9:12);                % State at time t3
        x3 = S3(1);                    % x-coordinate of position at time t3
        y3 = S3(2);                    % y-coordinate of position at time t3
        vx3 = S3(3);                   % x-coordinate of velocity at time t3
        vy3 = S3(4);                   % y-coordinate of velocity at time t3
        t3 = ti + (3-1)/(N-1)*(tf-ti); % Time t3


        % Propagate between requested times:
        [x_sol_Phi_12, PHIf_12, ~, ~, ~]  = propagation(ti, Si, t2, par, 'xyPBRFBP_STM_ROT');
        [x_sol_Phi_23, PHIf_23, ~, ~, ~]  = propagation(t2, S2, t3, par, 'xyPBRFBP_STM_ROT');
        [x_sol_Phi_34, PHIf_34, ~, ~, ~]  = propagation(t3, S3, tf, par, 'xyPBRFBP_STM_ROT');

        % Define the initial and final conditions:
        Psi1 = [(xi+mu)^2 + yi^2-ri^2; (xi+mu)*(vxi-yi)+yi*(vyi+xi+mu)];
        Psi2 = [(xf+mu-1)^2+yf^2-rf^2; (xf+mu-1)*(vxf-yf)+yf*(vyf+xf+mu-1)];

        % Continuity conditions:
        Concon1 = x_sol_Phi_12-S2;
        Concon2 = x_sol_Phi_23-S3;
        Concon3 = x_sol_Phi_34-Sf;

        % Define the equality constraints:
        ceq = [Psi1; Concon1; Concon2; Concon3; Psi2];

        % Define the inequality constraints:
        c = [E_Radii(1)^2/DU^2 - ((xi + mu)^2 + yi^2);     % initial distance from Earth
             M_Radii(1)^2/DU^2 - ((xi +mu - 1)^2 + yi^2);  % initial distance from Moon
             E_Radii(1)^2/DU^2 - ((x2 + mu)^2 + y2^2);     % intermediate distance from Earth
             M_Radii(1)^2/DU^2 - ((x2 +mu - 1)^2 + y2^2);  % intermediate distance from Moon
             E_Radii(1)^2/DU^2 - ((x3 + mu)^2 + y3^2);     % intermediate distance from Earth
             M_Radii(1)^2/DU^2 - ((x3 +mu - 1)^2 + y3^2);  % intermediate distance from Moon
             E_Radii(1)^2/DU^2 - ((xf + mu)^2 + yf^2);     % final distance from Earth
             M_Radii(1)^2/DU^2 - ((xf + mu - 1)^2 + yf^2); % final distance from Moon
             ti - tf];                                     % condition of tof>0
        

        if nargout > 2 % i.e. the gradients are requested
            
            % Define useful functions for the equality constraint gradient:
            dOM4dx = @(x,y,vx,vy,t) x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(t*ws))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(t*ws)))/(2*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2));
            dOM4dy = @(x,y,vx,vy,t) y - (ms*sin(t*ws))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(t*ws)))/(2*((x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);
            f = @(x, y, vx, vy, t) [vx; vy; dOM4dx(x,y,vx,vy,t) + 2*vy; dOM4dy(x,y,vx,vy,t) - 2*vx];
            
            Q_1_1 = -(N-1)/(N-1)*PHIf_12*f(xi, yi, vxi, vyi, ti) + (N-1-1)/(N-1)*f(x_sol_Phi_12(1), x_sol_Phi_12(2), x_sol_Phi_12(3), x_sol_Phi_12(4), t2);
            Q_2_1 = -(N-2)/(N-1)*PHIf_23*f(x2, y2, vx2, vy2, t2) + (N-2-1)/(N-1)*f(x_sol_Phi_23(1), x_sol_Phi_23(2), x_sol_Phi_23(3), x_sol_Phi_23(4), t3);
            Q_3_1 = -(N-3)/(N-1)*PHIf_34*f(x3, y3, vx3, vy3, t3) + (N-3-1)/(N-1)*f(x_sol_Phi_34(1), x_sol_Phi_34(2), x_sol_Phi_34(3), x_sol_Phi_34(4), tf);

            Q_1_N = -(1-1)/(N-1)*PHIf_12*f(xi, yi, vxi, vyi, ti) + 1/(N-1)*f(x_sol_Phi_12(1), x_sol_Phi_12(2), x_sol_Phi_12(3), x_sol_Phi_12(4), t2);
            Q_2_N = -(2-1)/(N-1)*PHIf_23*f(x2, y2, vx2, vy2, t2) + 2/(N-1)*f(x_sol_Phi_23(1), x_sol_Phi_23(2), x_sol_Phi_23(3), x_sol_Phi_23(4), t3);
            Q_3_N = -(3-1)/(N-1)*PHIf_34*f(x3, y3, vx3, vy3, t3) + 3/(N-1)*f(x_sol_Phi_34(1), x_sol_Phi_34(2), x_sol_Phi_34(3), x_sol_Phi_34(4), tf);
            
            R1 = [2*(xi+mu), 2*yi, 0, 0; vxi, vyi, xi+mu, yi];
            RN = [2*(xf+mu-1), 2*yf, 0, 0; vxf, vyf, xf+mu-1, yf];

            % Initialise the gradient of the equality constraint:
            gceq = zeros(16,18);

            % The gradient of the equality constraints is:
            gceq(1:2, 1:18) = [R1, zeros(2,4), zeros(2,4), zeros(2,4), zeros(2,1), zeros(2,1)];
            gceq(3:6, 1:18) = [PHIf_12, -eye(4), zeros(4,4), zeros(4,4), Q_1_1, Q_1_N];
            gceq(7:10, 1:18) = [zeros(4,4), PHIf_23, -eye(4), zeros(4,4), Q_2_1, Q_2_N];
            gceq(11:14, 1:18) = [zeros(4,4), zeros(4,4), PHIf_34, -eye(4), Q_3_1, Q_3_N];
            gceq(15:16, 1:18) = [zeros(2,4), zeros(2,4), zeros(2,4), RN, zeros(2,1), zeros(2,1)];
            gceq = gceq';

            % Define some useful quantities for the inequality constraint
            % gradient:
            S_1 = [-2*(xi+mu), -2*yi, 0, 0; -2*(xi+mu-1), -2*yi, 0, 0];
            S_2 = [-2*(x2+mu), -2*y2, 0, 0; -2*(x2+mu-1), -2*y2, 0, 0];
            S_3 = [-2*(x3+mu), -2*y3, 0, 0; -2*(x3+mu-1), -2*y3, 0, 0];
            S_N = [-2*(xf+mu), -2*yf, 0, 0; -2*(xf+mu-1), -2*yf, 0, 0];
            S_t = [1, -1];

            % Intialize gradient of the inequality constraints:
            gc = zeros(9, 18);

            % Compute the gradients of the inequality constraints:
            gc(1:2, 1:4) = S_1;
            gc(3:4, 5:8) = S_2;
            gc(5:6, 9:12) = S_3;
            gc(7:8, 13:16) = S_N;
            gc(9, 17:18) = S_t;
            gc = gc';

        end 
             
    else 
        error('Unknown method specified. Choose between ''simple'' and ''multiple_4''')
    end 
end

%--------------------------------------------------------------------------
function [orbitType, captureType] = characterise_orbits(S, par)
% DESCRIPTION:
% Characterizes the orbit and capture type based on the final values of angular
% momentum and Kepler energy.
% 
% PROTOTYPE:
% [orbitType, captureType] = characterise_orbits(S, par)
%
% INPUT:
% S[4x1]         State vector                           [-]
%   S(1)          x-coordinate of position             [DU]
%   S(2)          y-coordinate of position             [DU]
%   S(3)          x-velocity component                 [VU]
%   S(4)          y-velocity component                 [VU]
% par.mu[1x1]    Gravitational parameter of the central body [DU^3/TU^2]
%
% OUTPUT:
% orbitType[char]   Description of the orbit type, either 'direct capture', 
%                   'retrograde capture', or no residual rotational motion.
% captureType[char] Description of the capture type, either 'non-ballistically 
%                   captured', 'ballistically captured', or 'parabolic trajectory'.
    
    % Calculate final values
    h2 = (S(1) + par.mu-1)*(S(4)+S(1)+par.mu-1) - S(2)*(S(3)-S(2));       % Final value of angular momentum
    r_f = sqrt((S(1)+par.mu-1)^2 + S(2)^2);
    v_rel = sqrt((S(3)-S(2))^2 + (S(4)+S(1)+par.mu-1)^2);                 % Relavtive velocity 
    H_2 = v_rel^2/2 - par.mu/r_f;                                         % Final value of Kepler energy
    
    % Determine orbit type based on angular momentum
    if h2 > 0
        orbitType = sprintf('Since h2 = %f [DU*VU], the final orbit is a direct capture solution.', h2);
    elseif h2 < 0
        orbitType = sprintf('Since h2 = %f [DU*VU], the final orbit is a retrograde capture solution.', h2);
    else
        orbitType = sprintf('Since h2 = %f [DU*VU], there is no residual rotational motion relative to the primary body.', h2);
    end
    
    % Determine capture type based on Kepler energy
    if H_2 > 0
        captureType = sprintf('Since H2 = %f [VU^2], the SC is non-ballistically captured.', H_2);
    elseif H_2 < 0
        captureType = sprintf('Since H2 = %f [VU^2], the SC is ballistically captured.', H_2);
    else
        captureType = sprintf('Since H2 = %f [VU^2], the SC is on a parabolic trajectory.', H_2);
    end
end

%--------------------------------------------------------------------------
function [bodies] = nbody_init(labels)
% DESCRIPTION:
% Initializes planetary data for n-body propagation. Given a set of body labels,
% returns a cell array with structures containing the body labels and associated
% gravitational constants.
%
% PROTOTYPE:
% bodies = nbody_init(labels)
%
% INPUT:
% labels[1xn]      Cell array with labels of bodies (e.g., planets or barycenters)
%
% OUTPUT:
% bodies[1xn]      Cell array of structures with the following fields:
%   - bodies{i}.name -> Body label
%   - bodies{i}.GM   -> Gravitational constant [km^3/s^2]

    % Initialize output
    bodies = cell(size(labels));

    % Loop over labels
    for i = 1:length(labels)
        % Store body label
        bodies{i}.name = labels{i};
        % Store body gravitational constant
        bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
    end

end

%--------------------------------------------------------------------------
function [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
% DESCRIPTION
% Evaluates the right-hand-side of a Newtonian N-body propagator.
%
% PROTOTYPE
% [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
%
% INPUTS:
%   t      : [1,1] double - Ephemeris time (ET SPICE), seconds past J2000 (TDB).
%   x      : [6,1] double - Cartesian state vector (position and velocity) 
%                           with respect to the desired object.
%   bodies  : [1,n] cell-array - Cell array containing information about celestial 
%                           bodies, created with function nbody_init.
%   frame   : [1,x] string - Integration reference frame ('J2000' or 'ECLIPJ2000').
%   center  : [1,x] string - Name of the integration center (e.g., 'SSB').
%
% OUTPUTS:
%   dxdt   : [6,1] double - Right-hand side (RHS) vector, containing
%                           the position derivatives and gravitational accelerations.


    % Validate the integration frame
    if not(strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000'))
        error('Invalid integration reference frame. Select either J2000 or ECLIPJ2000.');
    end

    % Validate the integration center
    if not(any(strcmp(center, {cell2mat(bodies).name})))
        error('Invalid center selected. Select one of the bodies.');
    end

    % Initialize the right-hand side vector
    dxdt = zeros(6, 1);

    % Set the position derivative as the object's velocity
    dxdt(1:3) = x(4:6);

    % Extract the object's position from the state vector
    rr_b0_obj = x(1:3);

    % Extract the gravitational parameter (GM) of the central body
    GM0 = bodies{strcmp(center, {cell2mat(bodies).name})}.GM;

    % Compute the gravitational acceleration due to the central body
    dxdt(4:6) = -GM0 * rr_b0_obj / norm(rr_b0_obj)^3;

    % Loop over all bodies (except the central body)
    for i = 1:length(bodies)
        if strcmp(bodies{i}.name, center)
            continue; % Skip the central body
        end
        
        % Retrieve position and velocity of the i-th celestial body in the inertial frame
        rv_b0_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', center);
        
        % Calculate the position of the object with respect to the i-th celestial body
        rr_body_obj = rr_b0_obj - rv_b0_body(1:3);
        
        % Compute necessary variables for the non-inertial terms
        r = rr_b0_obj;             % Position vector of the object
        d = rr_body_obj;           % Position vector of the object with respect to the i-th body
        rho = rv_b0_body(1:3);     % Position vector of the i-th body in the inertial frame
        
        % Calculate q and f based on the positions
        q = dot(r, r - 2 * rho) / dot(rho, rho);
        f = q * (3 + 3 * q + q^2) / (1 + (1 + q)^(1.5));
        
        % Compute the gravitational contribution from the i-th body
        a_grav = -bodies{i}.GM * (1 / norm(d)^3) * (r + rho * f);
        
        % Sum up the acceleration contributions to the right-hand side vector
        dxdt(4:6) = dxdt(4:6) + a_grav;
    end
end

%--------------------------------------------------------------------------
function [penalty, theta] = thetafind(t, theta_i)
% DESCRIPTION:
% Calculates the angle between the Sun and the Moon's orbital plane, and computes
% the penalty based on the difference between the calculated and input angles.
%
% PROTOTYPE:
% [penalty, theta] = thetafind(t, theta_i)
%
% INPUT:
% t[1x1]          Time for which the angle is to be computed (Julian Date) [days]
% theta_i[1x1]     Initial angle value (rad) for penalty calculation
%
% OUTPUT:
% penalty[1x1]     Difference (penalty) between the computed and input angles [rad]
% theta[1x1]       Angle between the Sun and the Moon's orbital plane [rad]

    
        % Retrieve the ephemerides data of Moon and Sun:
        rr_sun = cspice_spkpos('Sun', t, 'ECLIPJ2000', 'NONE', 'EMB');       % No aberration correction
        rv_moon = cspice_spkezr('Moon', t, 'ECLIPJ2000', 'NONE', 'EMB');     % No aberration correction
        
        % Define the unitary vectors identifying the rotating frame:
        z_rot = cross(rv_moon(1:3), rv_moon(4:6)) / norm(cross(rv_moon(1:3), rv_moon(4:6)));
        x_rot = rv_moon(1:3)/norm(rv_moon);
        y_rot = cross(x_rot,-z_rot)/ norm(cross(x_rot,-z_rot));

        % Rotation matrix from ECLIP2000 to Moon orbit:
        T_ECLIP2000_2_Moon = [x_rot'; y_rot'; z_rot'];
 
        % Sun position vector rotated into the Moon orbital plane:
        rr_sun_p2 = T_ECLIP2000_2_Moon*rr_sun;
        
        % Find the desired angle:
        costheta = dot(rr_sun_p2, [1, 0, 0]) / (norm(rr_sun_p2));
        sintheta = dot(rr_sun_p2, [0, 1, 0]) / (norm(rr_sun_p2));
    
        theta = wrapTo2Pi(atan2(sintheta, costheta));

    % Calculate the penalty:
    penalty = angdiff(theta_i, theta);

end 

