% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #1, Exercise 1
% Author: Emanuele Gallo

% Prepare the workspace:
clc; clearvars; close all                                          % Clear command window, variables, and close figures

% Impose long format to better visualize results:
format long g                                                      % Use long format for better numerical precision

%% Exercise 1.1

% Define the (adimensional) gravitational constant:
mu = 0.012150;                                                    % Earth-Moon mass ratio [-]

% Define the (adimensional) position vectors:
r1_vect = @(x,y,z) [x+mu; y; z];                                  % Vector from the primary (Earth)
r1 = @(x,y,z) norm(r1_vect(x,y,z),2);                             % Distance from the primary (Earth) [-]
r2_vect = @(x,y,z) [x+mu-1; y; z];                                % Vector from the secondary (Moon)
r2 = @(x,y,z) norm(r2_vect(x,y,z),2);                             % Distance from the secondary (Moon) [-]

% Define the potential and its derivatives:
OM = @(x,y,z) (x^2+y^2)/2 + (1-mu)/r1(x,y,z) + mu/r2(x,y,z) + ... % Effective potential function
              mu*(1-mu)/2;                                        % Additional constant term
dOM_dx = @(x,y,z) x - (mu+x)*(1-mu)/(r1(x,y,z))^3 - ...           % x-derivative of the potential
                  mu*(x+mu-1)/(r2(x,y,z)^3);                     
dOM_dy = @(x,y,z) y - y*(1-mu)/(r1(x,y,z))^3 - ...                % y-derivative of the potential
                  mu*y/(r2(x,y,z)^3);                            
dOM_dz = @(x,y,z) -z*(1-mu)/(r1(x,y,z)^3) - ...                   % z-derivative of the potential
                  mu*z/(r2(x,y,z)^3);                            

% Define the Jacobi constant function:
J = @(S) 2*OM(S(1), S(2), S(3)) - (S(4)^2+S(5)^2+S(6)^2);         % Jacobi constant [-]

% Define intervals for L1, L2, and L3 based on physical reasoning: 
options = optimset('Display', 'off', 'TolX', 1e-15, 'TolFun', 1e-15); % Solver options for high precision

% L1: Between the Earth and Moon
L1_guess = [0.2, 1-2*mu];                                         % Initial guess for L1
L1_x = fzero(@(x) dOM_dx(x,0,0), L1_guess, options);              % Solve for L1's x-coordinate
L1 = [L1_x, 0, 0]';                                               % Define L1 position vector

% L2: Beyond the Moon
L2_guess = [1+mu, 2];                                             % Initial guess for L2
L2_x = fzero(@(x) dOM_dx(x,0,0), L2_guess, options);              % Solve for L2's x-coordinate
L2 = [L2_x, 0, 0]';                                               % Define L2 position vector

% L3: Opposite side of the Earth
L3_guess = [-2, -1+mu];                                           % Initial guess for L3
L3_x = fzero(@(x) dOM_dx(x,0,0), L3_guess, options);              % Solve for L3's x-coordinate
L3 = [L3_x, 0, 0]';                                               % Define L3 position vector

% L4 and L5: Triangular points
L45_x_guess = [0, 0.5];                                           % Initial guess for L4 and L5
L4_x = fzero(@(x) dOM_dx(x, sqrt(3)/2, 0), L45_x_guess, options); % Solve for L4's x-coordinate
L5_x = fzero(@(x) dOM_dx(x, -sqrt(3)/2, 0), L45_x_guess, options);% Solve for L5's x-coordinate
L4 = [L4_x, sqrt(3)/2, 0]';                                      % Define L4 position vector
L5 = [L5_x, -sqrt(3)/2, 0]';                                     % Define L5 position vector

% Calculate the Jacobi constant of all Lagrangian points:
C_L1 = J([L1; 0; 0; 0]);                                         % Jacobi constant at L1
C_L2 = J([L2; 0; 0; 0]);                                         % Jacobi constant at L2
C_L3 = J([L3; 0; 0; 0]);                                         % Jacobi constant at L3
C_L4 = J([L4; 0; 0; 0]);                                         % Jacobi constant at L4
C_L5 = J([L5; 0; 0; 0]);                                         % Jacobi constant at L5

% Display results:
fprintf('----------------Lagrangian points and related jacobi constant values--------------------------\n')
fprintf('L1: [%.10f, %.10f, %.10f]\n', L1);                       % Display L1 coordinates
fprintf('C_L1: %.10f\n', C_L1);                                   % Display L1 Jacobi constant
fprintf('L2: [%.10f, %.10f, %.10f]\n', L2);                       % Display L2 coordinates
fprintf('C_L2: %.10f\n', C_L2);                                   % Display L2 Jacobi constant
fprintf('L3: [%.10f, %.10f, %.10f]\n', L3);                       % Display L3 coordinates
fprintf('C_L3: %.10f\n', C_L3);                                   % Display L3 Jacobi constant
fprintf('L4: [%.10f, %.10f, %.10f]\n', L4);                       % Display L4 coordinates
fprintf('C_L4: %.10f\n', C_L4);                                   % Display L4 Jacobi constant
fprintf('L5: [%.10f, %.10f, %.10f]\n', L5);                       % Display L5 coordinates
fprintf('C_L5: %.10f\n ', C_L5);                            % Display L5 Jacobi constant
fprintf('Values expressed in Earth-Moon rotating frame centred @EMB\n\n\n');


% Plot:
x_plot = [linspace(-2, -mu, 300), linspace(-mu, 1-mu, 300), linspace(1-mu, 2, 300)];
dOM_plot1 = zeros(1, length(x_plot));
OM_plot1 = zeros(1, length(x_plot));
dOM_plot2 = zeros(1, length(x_plot));
OM_plot2 = zeros(1, length(x_plot));
dOM_plot3 = zeros(1, length(x_plot));
OM_plot3 = zeros(1, length(x_plot));

for i=1:length(x_plot)
    dOM_plot1(i) = dOM_dx(x_plot(i), 0, 0);
    OM_plot1(i) = OM(x_plot(i), 0, 0);
    dOM_plot2(i) = dOM_dx(x_plot(i), sqrt(3)/2, 0);
    OM_plot2(i) = OM(x_plot(i), sqrt(3)/2, 0);
    dOM_plot3(i) = dOM_dx(x_plot(i), -sqrt(3)/2, 0);
    OM_plot3(i) = OM(x_plot(i), -sqrt(3)/2, 0);

    if dOM_plot1(i)>40
        dOM_plot1(i) = NaN;
    elseif dOM_plot2(i)>20
        dOM_plot2(i) = Nan;
    elseif dOM_plot3(i)>20
        dOM_plot3(i) = NaN;
    end
end

figure; hold on; grid on;

% Plot Lagrange points
plot(L1(1), L1(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Collinear points');
plot(L2(1), L2(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility','off');
plot(L3(1), L3(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility','off');
plot(L4(1), L4(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Triangular points');
plot(L5(1), L5(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility','off');

% Create Earth and Moon representations
theta = linspace(0, 2*pi, 100); % For circles

% Earth
earth_radius = 0.07;
earth_x = earth_radius * cos(theta) - mu; % Earth position
earth_y = earth_radius * sin(theta);
fill(earth_x, earth_y, [0.2 0.4 1], 'EdgeColor', 'none'); % Blue Earth
text(-mu, -0.12, 'Earth', 'Interpreter','latex', 'FontSize', 8, 'HorizontalAlignment', 'center');

% Moon
moon_radius = 0.05;
moon_x = moon_radius * cos(theta) + (1-mu); % Moon position
moon_y = moon_radius * sin(theta);
fill(moon_x, moon_y, [0.6 0.6 0.6], 'EdgeColor', 'none'); % Grey Moon
text(1-mu, -0.12, 'Moon', 'Interpreter','latex', 'FontSize', 8, 'HorizontalAlignment', 'center');

% Add circles around Earth and Moon
r = 1; % Radius of influence zones
circle_earth_x = r*cos(theta) - mu; % Circle around Earth
circle_earth_y = r*sin(theta);
circle_moon_x = r*cos(theta) + (1-mu); % Circle around Moon
circle_moon_y = r*sin(theta);
plot(circle_earth_x, circle_earth_y, 'b--', 'LineWidth', 1, 'DisplayName', 'Earth Influence Zone');
plot(circle_moon_x, circle_moon_y, 'k--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Moon Influence Zone');

% Add text labels for the circles
text(-mu - 0.2, -r - 0.1, 'Circle with $R=1$  centered on Earth', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center');
text(1-mu + 0.2, r + 0.1, 'Circle with $R=1$  centered on Moon', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center');

% Add labels for Lagrange points
text(L1(1)-0.06, L1(2) + 0.1, '$L_1$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'left');
text(L2(1)-0.06, L2(2) + 0.1, '$L_2$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'left');
text(L3(1) + 0.06, L3(2), '$L_3$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'left');
text(L4(1)-0.06, L4(2) + 0.1, '$L_4$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'left');
text(L5(1)-0.06, L5(2) + 0.1, '$L_5$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'left');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');

% Title
title('Lagrangian points and Earth and Moon @EMB Earth-Moon rotating frame', 'FontSize', 10);

% Set axis limits and grid
axis([-1.5 2.5 -2 2]);
grid on;
hold off;

% Plot potential
figure; hold on; grid on;

% Plot potentials
plot(x_plot, OM_plot1, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Potential $\Omega(x, 0, 0)$');
plot(x_plot, OM_plot2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Potential $\Omega(x, \pm \frac{\sqrt{3}}{2}, 0)$');

% Add labels
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('$\Omega$ [-]', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');

% Add vertical lines and text annotations
xline(-mu, 'LineStyle', '--', 'LineWidth', 1, 'Color', [0, 0, 0], 'DisplayName', 'Primary bodies');
text(-mu-0.2, 0.5, 'Earth', 'FontSize', 10, 'Color', [0 0 0], 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
xline(1-mu, 'LineStyle', '--', 'LineWidth', 1, 'Color', [0, 0, 0], 'HandleVisibility', 'off');
text(1-mu, 0.5, 'Moon', 'FontSize', 10, 'Color', [0 0 0], 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Add Lagrange points
xline(L1(1), 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0 0], 'DisplayName', 'Collinear points');
xline(L2(1), 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0 0], 'HandleVisibility', 'off');
xline(L3(1), 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0 0], 'HandleVisibility', 'off');
xline(L4(1), 'LineStyle', '--', 'LineWidth', 1, 'Color', [0 0 1], 'DisplayName', 'Triangular points');
text(L1(1) - 0.1, 0, '$L_1$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', [1 0 0], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(L2(1) + 0.1, 0, '$L_2$', 'Interpreter', 'latex', 'Color', [1 0 0], 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(L3(1) - 0.1, 0, '$L_3$', 'Interpreter', 'latex', 'Color', [1 0 0], 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(L4(1), -0.3, '$L_4 \equiv L_5$', 'Interpreter', 'latex', 'Color', [0 0 1], 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Add legend
legend('show', 'Interpreter', 'latex', 'FontSize', 9, 'location', 'northwest');
title('Potential', 'Interpreter', 'latex', 'FontSize', 10)

% Set axis limits
axis([-2 2 -0.5 5]);

hold off;




%% Exercise 1.2

% Initial guess:
x0 = 1.068792441776;                  % Initial x-coordinate
y0 = 0;                               % Initial y-coordinate
z0 = 0.071093328515;                  % Initial z-coordinate
vx0 = 0;                              % Initial velocity in x-direction
vy0 = 0.319422926485;                 % Initial velocity in y-direction
vz0 = 0;                              % Initial velocity in z-direction
S0 = [x0, y0, z0, vx0, vy0, vz0]';    % Initial state vector

% Propagation time set to 1 orbit from null initial time:
t0 = 0;                               % Initial time
T_orbit = 2*pi;                       % Orbital period (1 full orbit)
tf = T_orbit;                         % Final time of propagation

% Reference value of the Jacobi constant:
C_ref = 3.09;                         % Target Jacobi constant value

% Numerical aid values:
Nmax = 100;                           % Maximum iterations to avoid divergence
tol = [1e-6; 1e-9; 1e-9; 1e-7];       % Convergence tolerances

% Apply the variational approach of STM to correct the initial conditions:
[S0_new, te, corr, ~] = diff_corr_STM(S0, C_ref, 0, tf, mu, Nmax, tol);

fprintf('------------------------------Differential correction--------------------\n')
fprintf('The corrected initial state for C = %.3f is S0_corrected =\n', C_ref);
fprintf('[%.10f] \n ', S0_new');      % Display corrected initial state
fprintf('Values expressed in Earth-Moon rotating frame centred @EMB\n\n\n');

% Check over the correctness of the correction:
if corr == 1                          % If correction was successful
    % Propagate the orbit with corrected initial conditions:
    [xf_fb, PHIf_f, tf_f, xx_f, tt_f] = propagation(0, S0_new, te, mu, true); % Forward propagation
    [xf_b, PHIf_b, tf_b, xx_b, tt_b] = propagation(0, S0_new, -te, mu, true); % Backward propagation
    
    % Plot the propagated orbit with corrected initial conditions:
    figure();
    plot3(xx_f(:,1), xx_f(:,2), xx_f(:,3), 'LineWidth', 2, 'Color', [0, 0, 1]); % Forward trajectory (blue) 
    hold on;                            
    plot3(xx_b(:,1), xx_b(:,2), xx_b(:,3), 'LineWidth', 2, 'Color', [0, 0, 1]); % Backward trajectory (blue) 

    % Plot central bodies:
    plot3(1-mu, 0, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'MarkerEdgeColor', 'b'); % Moon's position
    text(1-mu + 0.008, 0.008, 0.008, 'Moon', 'Interpreter', 'latex', 'FontSize', 12, 'HorizontalAlignment', 'center'); % Label Moon

    % Plot L2 point:
    plot3(L2(1), L2(2), L2(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); % L2 point
    text(L2(1) + 0.008, L2(2) + 0.008, L2(3) + 0.008, '$L_{2}$', 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'center'); % Label L2
    
    % Add 'Orbit' label:
    text(1.12927775217553, 0.0814678087635405, 0.0145123316629201, 'Orbit', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Label along trajectory

    % Plot formatting:
    xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 50, 'FontWeight', 'bold'); % X-axis label
    ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 50, 'FontWeight', 'bold'); % Y-axis label
    zlabel('$z$ [-]', 'Interpreter', 'latex', 'FontSize', 50, 'FontWeight', 'bold'); % Z-axis label
    title('Orbit with Corrected Initial Conditions (@EMB Earth-Moon rotating frame)', 'Interpreter', 'latex', 'FontSize', 2, 'FontWeight', 'bold'); % Plot title
    
    % Additional settings:
    grid on;                          % Enable grid
    axis equal;                       % Equal scaling for axes
    view(3);                          % Set 3D view perspective
    set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Style axes

    % Background and lighting:
    set(gca, 'Color', [1 1 1]);       % Set background color
    camlight left;                    % Add lighting from the left
    lighting phong;                   % Use Phong lighting for smoothness
    hold off;                         % End plotting
end


%% Exercise 1.3

% Useful jacobi constant values:
C_f = 3.04;                          % Final Jacobi constant value wanted
C_vect = linspace(C_ref, C_f, 20);   % Vector of Jacobi constants from reference to final value

% Propagation time set to 1 orbit:
T_orbit = 2*pi;                      % Period of one orbit
tf = T_orbit;                        % Final time for propagation

% Initialize the variables:
S_vect = zeros(6, 1, length(C_vect));        % Preallocate array for state vectors
S_vect(1:6, 1, 1) = S0_new;                  % Assign corrected initial state to first element

% Create new figure to avoid superposition with the previous plots:
figure()
set(gcf, 'Color', 'w');                      % Set white background for the figure

% Colormap setup:
num_colors = length(C_vect);                 % Number of orbits to color
colors = parula(num_colors);                 % Assign colormap (parula)

% Legend handles array:
legend_handles = [];                         % Initialize array for legend entries

% Plot the L2 point (no legend entry for this):
plot3(L2(1), L2(2), L2(3), 'r.', 'MarkerSize', 20);               % L2 point as a red dot
hold on
text(L2(1), L2(2) + 0.01, L2(3) + 0.01, '$L_{2}$', ...            % Label the L2 point
    'Interpreter', 'latex', 'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'center');

% Plot moon (no red circle, only black marker):
plot3(1-mu, 0, 0, '.', 'MarkerSize', 20, 'Color', 'k');           % Moon as a black dot
text(1-mu, 0.01, 0.01, 'Moon', 'Interpreter', 'latex', ...        % Label the Moon
    'FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'center');

% Initial orbit:
h_initial = plot3(xx_f(:,1), xx_f(:,2), xx_f(:,3), ...            % Forward propagation (initial orbit)
    'LineWidth', 3, 'DisplayName', sprintf('Initial Orbit with C = %.3f', C_vect(1)), ...
    'Color', colors(1, :));
legend_handles = [legend_handles, h_initial];                    % Add handle to legend
plot3(xx_b(:,1), xx_b(:,2), xx_b(:,3), 'LineWidth', 3, ...        % Backward propagation (initial orbit)
    'Color', colors(1, :));                                       % Same color as forward

% Loop for intermediate orbits:
intermediate_orbit_idx = 4;                                       % Plot every 4th orbit
for i = 2:length(C_vect)
    % Correct initial conditions:
    [S_vect(1:6,1,i), te(i), corr, xx_f3] = diff_corr_STM(S_vect(1:6,1,i-1), C_vect(i), t0, tf, mu, Nmax, tol);

    % Check the correction and plot if successful:
    if corr == 1
        if mod(i, intermediate_orbit_idx) == 0 && i ~= length(C_vect) % Plot only selected orbits
            [~, ~, tf_b, xx_b3, tt_b] = propagation(0, S_vect(1:6,1,i), -te(i), mu, true);

            % Assign color for the current orbit:
            orbit_color = colors(i, :);

            % Plot forward propagation:
            h_intermediate = plot3(xx_f3(:,1), xx_f3(:,2), xx_f3(:,3), 'LineWidth', 3, ...
                'Color', orbit_color, 'DisplayName', sprintf('Intermediate Orbit with C = %.3f', C_vect(i)));
            legend_handles = [legend_handles, h_intermediate];    % Add handle to legend

            % Plot backward propagation (no legend entry):
            plot3(xx_b3(:,1), xx_b3(:,2), xx_b3(:,3), 'LineWidth', 3, 'Color', orbit_color);
        elseif i == length(C_vect)                                % Final orbit
            [~, ~, tf_b, xx_b3, tt_b] = propagation(0, S_vect(1:6,1,i), -te(i), mu, true);
        end
    end
end

% Final orbit (ensure xx_f3 is valid):
if exist('xx_f3', 'var') && ~isempty(xx_f3)
    final_color = colors(end, :);                                % Assign final color
    h_final = plot3(xx_f3(:, 1), xx_f3(:, 2), xx_f3(:, 3), ...   % Plot forward final orbit
        'LineWidth', 3, 'Color', final_color, ...
        'DisplayName', sprintf('Final Orbit with C = %.3f', C_vect(end)));
    legend_handles = [legend_handles, h_final];                  % Add handle to legend
    plot3(xx_b3(:, 1), xx_b3(:, 2), xx_b3(:, 3), 'LineWidth', 3, 'Color', final_color); % Backward final orbit
else
    warning('Final orbit data (xx_f3) is empty or invalid. Skipping final plot.'); % Warning for invalid data
end

% Show legend with selected handles only:
legend(legend_handles, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10);

% Grid, axis labels, and scaling:
grid on
axis equal
xlabel('$x[-]$', 'Interpreter', 'latex', 'FontSize', 20);         % x-axis label
ylabel('$y[-]$', 'Interpreter', 'latex', 'FontSize', 20);         % y-axis label
zlabel('$z[-]$', 'Interpreter', 'latex', 'FontSize', 20);         % z-axis label
title(sprintf('Orbit with Corrected Initial Conditions for varying C (Earth-Moon rotating frame @EMB)'), 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold'); % Plot title
hold off                                                          % End plotting

% Display corrected final state:
fprintf('----------------------Numerical continuation-----------------------\n')
fprintf('The corrected initial state for C = %.3f is S0_corrected = \n', C_vect(end));
fprintf('[%.10f] \n', S_vect(1:6,1,end)');                        % Print final state vector
fprintf('Values expressed in Earth-Moon rotating frame centred @EMB\n\n\n');



%% Functions

% -------------------------------------------------------------------------
function [dSdt] = xyzCR3BP_STM(~, S, mu)
% DESCRIPTION:
% Computes the Right-Hand Side (RHS) of the Circular Restricted 3-Body 
% Problem (CR3BP) equations, including the State Transition Matrix (STM).
% 
% PROTOTYPE:
% [dSdt] = xyzCR3BP_STM(~, S, mu)
%
% INPUT:
% ~       Placeholder for time (not used)                     [-]
% S[42x1] State vector, including position, velocity, 
%         and STM elements                                    [-]
% mu[1x1] Gravitational parameter ratio of the two bodies     [-]
%
% OUTPUT:
% dSdt[42x1] Time derivative of the state vector              [-]

    % Extract variables from the state vector:
    x = S(1);       % x-coordinate
    y = S(2);       % y-coordinate
    z = S(3);       % z-coordinate
    vx = S(4);      % x-velocity
    vy = S(5);      % y-velocity
    vz = S(6);      % z-velocity

    Phi = reshape(S(7:end), 6, 6);                           % STM (6x6)

    % Compute distances to the primary bodies:
    r1 = sqrt((x + mu)^2 + y^2 + z^2);                       % Distance to body 1
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);                   % Distance to body 2

    % Compute potential derivatives:
    dOMdx = x - (1-mu)/r1^3*(mu + x) + mu/r2^3*(1-mu - x);   % Partial derivative w.r.t. x
    dOMdy = y - (1-mu)/r1^3*y - mu/r2^3*y;                  % Partial derivative w.r.t. y
    dOMdz = - (1-mu)/r1^3*z - mu/r2^3*z;                    % Partial derivative w.r.t. z


    % Assemble the matrix A = dfdx
    dfdx = [0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  1, 0, 0;
            0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 1, 0;
            0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 0, 1;
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1,                                                                     (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),                                                                 (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),  0, 2, 0;
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1,                                                                                   (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), -2, 0, 0;
            (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),                                                                                       (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2),  0, 0, 0];
 
 
    % Compute derivative of the STM:
    Phidot = dfdx * Phi;                                    % STM derivative

    % Assemble the RHS of the system:
    dSdt = zeros(42, 1);                                     % Initialize derivative vector
    dSdt(1:3) = [vx; vy; vz];                                % Position derivatives
    dSdt(4) = dOMdx + 2*vy;                                  % x-velocity derivative
    dSdt(5) = dOMdy - 2*vx;                                  % y-velocity derivative
    dSdt(6) = dOMdz;                                         % z-velocity derivative
    dSdt(7:end) = Phidot(:);                                 % Flattened STM derivative

end

% -------------------------------------------------------------------------
function [xf, PHIf, tf, xx, tt] = propagation(t0, x0, tf, mu, varargin)
% DESCRIPTION:
% Propagates the state vector and State Transition Matrix (STM) in the 
% Circular Restricted 3-Body Problem (CR3BP) dynamics until a specified 
% event or final time is reached.
% 
% PROTOTYPE:
% [xf, PHIf, tf, xx, tt] = propagation(t0, x0, tf, mu, varargin)
%
% INPUT:
% t0[1x1]       Initial time                                [s]
% x0[6x1]       Initial state vector (position and velocity) [-]
% tf[1x1]       Final time                                  [s]
% mu[1x1]       Gravitational parameter ratio of the two 
%               bodies                                      [-]
% varargin      Optional input for event handling
%               - evtFlag[1x1]: Enables/disables y=0 
%                 plane crossing detection                  [logical]
% 
% OUTPUT:
% xf[6x1]       Final state vector                         [-]
% PHIf[6x6]     Final State Transition Matrix (STM)        [-]
% tf[1x1]       Final time reached                         [s]
% xx[Nx6]       State vectors at all integration steps     [-]
% tt[Nx1]       Time vector for the propagation            [s]

    if nargin > 4                                             % Check for optional input
        evtFlag = varargin{1};                               % Event handling flag
    else
        evtFlag = true;                                      % Default: enable events
    end

    % Preintegration operations:
    Phi0 = eye(6);                                           % Initialize STM as identity matrix
    x0Phi0 = [x0; Phi0(:)];                                  % Concatenate state and STM

    % Set integration options, including tolerances and event function:
    options_STM = odeset('reltol', 1e-12, ...
                         'abstol', [1e-10*ones(1,3), ...     % Position tolerance
                                   1e-13*ones(1,3), ...      % Velocity tolerance
                                   1e-13*ones(1,36)], ...    % STM tolerance
                         'Events', @(x, y) xz_plane_crossing(x, y, evtFlag));

    % Integrate the dynamics using ode78:
    [tt, xx] = ode78(@(t, x) xyzCR3BP_STM(t, x, mu), [t0 tf], x0Phi0, options_STM);

    % Extract the final state and STM:
    xf = xx(end, 1:6)';                                      % Final state vector
    PHIf = reshape(xx(end, 7:end), 6, 6);                    % Final STM
    tf = tt(end);                                            % Final time
end

%--------------------------------------------------------------------------
function [value, isterminal, direction] = xz_plane_crossing(~, y, isTerminal)
% DESCRIPTION:
% Defines the conditions for detecting the crossing of the xz-plane during
% numerical integration.
%
% PROTOTYPE:
% [value, isterminal, direction] = xz_plane_crossing(~, y, isTerminal)
%
% INPUT:
% ~             Placeholder for time variable (unused)                  [-]
% y[6x1]        Current state vector (position and velocity)            [-]
% isTerminal[1x1] Flag indicating whether the integration should stop   [-]
%
% OUTPUT:
% value[1x1]    Value indicating the condition for plane crossing       [-]
%               (0 when the y-component is zero).
% isterminal[1x1] Flag to stop integration when crossing condition is met [-]
% direction[1x1] Direction of crossing (0 for any direction)            [-]

    value = y(2);          % Check y-coordinate to determine crossing
    isterminal = isTerminal; % Use provided terminal condition
    direction = 0;          % Allow detection in both directions
end

% -----------------------------------------------------------------------
function [S0_new, te, corr, xx] = diff_corr_STM(S0, C_ref, t0, tf, mu, Nmax, tol)
% DESCRIPTION:
% Implements differential correction using State Transition Matrix (STM)
% to adjust initial conditions for achieving a periodic halo orbit with
% a specified Jacobi constant.
%
% PROTOTYPE:
% [S0_new, te, corr, xx] = diff_corr_STM(S0, C_ref, t0, tf, mu, Nmax, tol)
%
% INPUT:
% S0[6x1]       Initial state vector                                   [-]
% C_ref[1x1]    Target Jacobi constant                                 [-]
% t0[1x1]       Initial time                                           [-]
% tf[1x1]       Final propagation time                                [-]
% mu[1x1]       Gravitational parameter of the CR3BP system           [-]
% Nmax[1x1]     Maximum number of iterations for convergence (default: 100) [-]
% tol[4x1]      Tolerance vector: [tol_y, tol_vx, tol_vz, tol_C]      [-]
%
% OUTPUT:
% S0_new[6x1]   Corrected initial state                               [-]
% te[1x1]       Event time for y=0 crossing                           [-]
% corr[1x1]     Convergence status (1: successful, 0: not needed, -1: failed) [-]
% xx[??x6]      State evolution over time for plotting purposes       [-]

    % Impose default values if not provided
    if nargin < 4
        Nmax = 100;                        % Default max iterations
        tol = [1e-6; 1e-9; 1e-9; 1e-7];    % Default tolerances
    end

    % Extract initial state components
    x0 = S0(1); y0 = S0(2); z0 = S0(3);    % Initial positions
    vx0 = S0(4); vy0 = S0(5); vz0 = S0(6);% Initial velocities

    % Compute initial Jacobi constant (C0)
    r1_0 = norm([x0 + mu; y0; z0]);        % Distance to primary
    r2_0 = norm([x0 + mu - 1; y0; z0]);   % Distance to secondary
    OM_0 = (x0^2 + y0^2) / 2 + (1 - mu) / r1_0 + mu / r2_0 + mu * (1 - mu) / 2;
    C0 = 2 * OM_0 - (vx0^2 + vy0^2 + vz0^2);

    % Check if correction is necessary
    if abs(C0 - C_ref) <= 1e-15            % If error is negligible
        fprintf('For C = %f, correction is not necessary.\n', C_ref);
        S0_new = S0; corr = 0;             % No correction needed
        return;
    end

    % Initialize correction variables
    err_y = 1; err_vxf = 1; err_vzf = 1; err_C = 1; % Errors
    iter = 0;                                   % Iteration counter
    x0_new = x0; z0_new = z0; vy0_new = vy0; tf_new = tf; % State updates

    % Set propagation targets
    vxf_ref = 0; vzf_ref = 0; y_ref = 0;

    % Iterative correction loop
    while (abs(err_y) > tol(1) || abs(err_vxf) > tol(2) || ...
           abs(err_vzf) > tol(3) || abs(err_C) > tol(4)) && iter < Nmax

        iter = iter + 1; % Increment iteration count

        % Propagate up to event time
        [xf, PHIf, te, xx, ~] = propagation(t0, [x0_new, y0, z0_new, vx0, vy0_new, vz0]', tf_new, mu, true);

        % Compute new Jacobi constant and linear dynamics
        r1_C = norm([x0_new + mu; y0; z0_new]);
        r2_C = norm([x0_new + mu - 1; y0; z0_new]);
        dOM_dx_C = x0_new - (mu + x0_new) * (1 - mu) / r1_C^3 - mu * (x0_new + mu - 1) / r2_C^3;
        dOM_dy_C = y0 - y0 * (1 - mu) / r1_C^3 - mu * y0 / r2_C^3;
        dOM_dz_C = -z0_new * (1 - mu) / r1_C^3 - mu * z0_new / r2_C^3;
        C_dyn_value = [2 * dOM_dx_C, 2 * dOM_dy_C, 2 * dOM_dz_C, -2 * vx0, -2 * vy0_new, -2 * vz0];
        OM_C = (x0_new^2 + y0^2) / 2 + (1 - mu) / r1_C + mu / r2_C + mu * (1 - mu) / 2;
        C = 2 * OM_C - (vx0^2 + vy0_new^2 + vz0^2);

        % Define useful quantities
        r1_f = norm([xf(1) + mu; xf(2); xf(3)]);
        r2_f = norm([xf(1) + mu - 1; xf(2); xf(3)]);
        dOM_dx_f = xf(1) - (mu + xf(1)) * (1 - mu) / r1_f^3 - mu * (xf(1) + mu - 1) / r2_f^3;
        dOM_dz_f = -xf(3) * (1 - mu) / r1_f^3 - mu * xf(3) / r2_f^3;

        % Coefficient matrix
        B = [PHIf(2,1), PHIf(2,3), PHIf(2,5), xf(5);
            PHIf(4,1), PHIf(4,3), PHIf(4,5), dOM_dx_f + 2 * xf(5);
            PHIf(6,1), PHIf(6,3), PHIf(6,5), dOM_dz_f;
            C_dyn_value(1), C_dyn_value(3), C_dyn_value(5), 0];

        % Known values
        b = [y_ref - xf(2); vxf_ref - xf(4); vzf_ref - xf(6); C_ref - C];

        % Solve to get the corrections:
        delta = B \ b;

        % Apply corrections
        x0_new = x0_new + delta(1);
        z0_new = z0_new + delta(2);
        vy0_new = vy0_new + delta(3);
        tf_new = tf_new + delta(4);

        % Update errors
        err_y = y_ref - xf(2);
        err_vxf = vxf_ref - xf(4);
        err_vzf = vzf_ref - xf(6);
        err_C = C_ref - C;
    end

    % Check convergence status
    S0_new = [x0_new; y0; z0_new; vx0; vy0_new; vz0];
    if iter < Nmax
        corr = 1; % Successful convergence
    else
        corr = -1; % Failed convergence
        fprintf('For C = %f, correction did not converge.\n', C_ref);
    end
end
