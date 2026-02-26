% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2, Exercise 1
% Author: Emanuele Gallo

% Prepare the workspace:
clear 
close all
clc

% Impose long format to better visualize results:
format long g

% Load the kernels:
cspice_furnsh('assignment02.tm')

%% Exercise 1.1

% Constants to solve the PCRTBP:
par.mu = 0.0121505842699404; % gravitational parameter
par.ms = 3.28900541e5;       % [-]: Scaled mass of the Sun
par.rho = 3.88811143e2;      % [-]: Scaled Sun-(Earth+Moon) distance
par.ws = -9.25195985e-1;     % [-]: Scaled angular velocity of the Sun
par.wem = 2.66186135e-1;     % [s^(-1)]: Earth-Moon angular velocity 
par.lem = 3.84405e8;         % [m]: Earth-Moon distance
par.hi = 167;                % [km]: Altitude of departure orbit
par.hf = 100;                % [km]: Altitude of arrival orbit
par.DU = 3.84405000e5;       % [km]: Distance Unit
par.TU = 4.34256461;         % [days]: Time Unit
par.VU = 1.02454018;         % [km/s]: Velocity Unit
par.T_orbit = 2*pi;          % [-]: Orbital period (adimensional)

% Extract the Moon radius from the kernels:
par.M_Radii = cspice_bodvrd('Moon', 'RADII',3);    % [km]: Radius of the Moon
R_M = par.M_Radii(1) / par.DU;                     % [DU]: Adimensionalised radius of the moon

% Take the values from the table:
ri = [-0.011965533749906, -0.017025663128129];     % [DU]:  Initial position mean state for satellite 
vi = [10.718855256727338, 0.116502348513671];      % [VU]:  Initial velocity mean state for satellite 
xi = [ri, vi];                                     % [DU, VU]: Initial state
ti = 1.282800225339865;                            % [TU]: Initial time 
tf = 9.595124551366348;                            % [TU]: Final time
t_vect = linspace(ti, tf, 5);                      % [TU]: Vector of times

% Initial covariance matrix
P0 = [1.041e-15, 6.026e-17, 5.647e-16, 4.577e-15; 
      6.026e-17, 4.287e-18, 4.312e-17, 1.855e-16;
      5.647e-16, 4.312e-17, 4.432e-16, 1.455e-15;
      4.577e-15, 1.855e-16, 1.455e-15, 2.822e-14];

% Unscented Transform parameters
alpha = 1;
beta = 2;

% Call LinCov method:
[E_LC, P_LC] = LinCov(xi', P0, t_vect, par);
mean_LC = E_LC(:, end);                         % Final mean (Linear Covariance)
cov_LC = P_LC(:,:,end);                         % Final covariance (Linear Covariance)

% Call Unscented Transform method:
[E_UT, P_UT] = UnTr(xi', P0, t_vect, alpha, beta, par);
mean_UT = E_UT(:,end);                                  % Final mean (Unscented Transform)
cov_UT = P_UT(:,:,end);                                 % Final covariance (Unscented Transform)

% Calculate the radius of the circular orbit (in DU):
r_orbit = (par.M_Radii(1) + par.hf) / par.DU; % [DU]: Radius of orbit

% Generate points for the circular orbit:
theta = linspace(0, 2*pi, 1000); % Angle values for the circle
x_orbit = r_orbit * cos(theta) + 1 - par.mu; % x-coordinates in rotating frame
y_orbit = r_orbit * sin(theta); % y-coordinates in rotating frame

% Plot:
options_plot_LC = struct('color', 'r', 'alpha', 0.2, 'linewidth', 1, 'linestyle', '-', 'DisplayName', 'LinCov 3$\sigma$ ellipse');
options_plot_UT = struct('color', 'b', 'alpha', 0.2, 'linewidth', 1, 'linestyle', '-', 'DisplayName', 'UT 3$\sigma$ ellipse');

% Plot the main figure
figure;
hold on;
grid on;

% Plot the reference circular orbit:
plot(x_orbit, y_orbit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference final orbit');

% Plot Linear Covariance mean and uncertainty ellipse:
plot(mean_LC(1), mean_LC(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'LinCov Mean');
plot_uncertainty_ellipse(mean_LC(1:2), cov_LC(1:2, 1:2), 3, options_plot_LC);

% Plot Unscented Transform mean and uncertainty ellipse:
plot(mean_UT(1), mean_UT(2), 'bx', 'MarkerFaceColor', 'b', 'MarkerSize', 12, 'DisplayName', 'UT Mean');
plot_uncertainty_ellipse(mean_UT(1:2), cov_UT(1:2, 1:2), 3, options_plot_UT);

xlabel('$x [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
title('General view of the comparison of Position Mean and Uncertainty Ellipses (@EMB Earth-Moon rotating at $t_f$)', 'Interpreter', 'latex', 'FontSize', 8);
axis equal;
hold off;

% Plot the main figure
figure;
hold on;
grid on;

% Plot the reference circular orbit:
plot(x_orbit, y_orbit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference final orbit');

% Plot Linear Covariance mean and uncertainty ellipse:
plot(mean_LC(1), mean_LC(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'LinCov Mean');
plot_uncertainty_ellipse(mean_LC(1:2), cov_LC(1:2, 1:2), 3, options_plot_LC);

% Plot Unscented Transform mean and uncertainty ellipse:
plot(mean_UT(1), mean_UT(2), 'bx', 'MarkerFaceColor', 'b', 'MarkerSize', 12, 'DisplayName', 'UT Mean');
plot_uncertainty_ellipse(mean_UT(1:2), cov_UT(1:2, 1:2), 3, options_plot_UT);

xlabel('$x [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
title('Zoomed-in view of the comparison of Position Mean and Uncertainty Ellipses (@EMB Earth-Moon rotating at $t_f$)', 'Interpreter', 'latex', 'FontSize', 8);
axis equal;
axis([mean_UT(1)-0.00165 mean_UT(1)+0.00165 mean_UT(2)-0.00165 mean_UT(2)+0.00165])
hold off;



%% Exercise 1.2

% Monte Carlo simulation:
n_samples = 5000;                                              % Decide the number of samples for the Monte Carlo analysis
[E_MC, P_MC, Sf_MC] = MC(xi', P0, t_vect, par, n_samples);     % Run the Monte Carlo analysis
mean_MC = E_MC(:, end);                                        % Final mean (Monte Carlo method)
cov_MC = P_MC(:,:,end);                                        % Final covariance matrix (Monte Carlo method)

% Results plot 1:
figure; hold on; grid on;
plot(Sf_MC(1, :), Sf_MC(2, :), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.7, 0.7, 0.7], ...
     'MarkerEdgeColor', 'none', 'DisplayName', 'MC Samples at final time', 'MarkerSize', 2);                           % Monte Carlo simulation points
plot(x_orbit, y_orbit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference final orbit');
options_plot_LC = struct('color', 'r', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '-', 'DisplayName', 'LinCov 3$\sigma$ ellipse');
plot(mean_LC(1), mean_LC(2), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'LinCov Mean'); % Mean of LinCov
plot_uncertainty_ellipse(mean_LC(1:2), cov_LC(1:2, 1:2), 3, options_plot_LC);                         % Ellipse of LinCov
options_plot_UT = struct('color', 'b', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '--', 'DisplayName', 'UT 3$\sigma$ ellipse');
options_plot_MC = struct('color', 'g', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '-', 'DisplayName', 'MC 3$\sigma$ ellipse');
plot(mean_MC(1), mean_MC(2), 'g+', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'MC Mean');     % Mean of MC
plot_uncertainty_ellipse(mean_MC(1:2), cov_MC(1:2, 1:2), 3, options_plot_MC);                         % Ellipse of MC
plot(mean_UT(1), mean_UT(2), 'bx', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'UT Mean');     % Mean of UT
plot_uncertainty_ellipse(mean_UT(1:2), cov_UT(1:2, 1:2), 3, options_plot_UT);                         % Ellipse of UT
xlabel('$x \, [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y \, [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
axis equal
title('General view of comparison of Position Mean and Uncertainty Ellipses (@EMB Earth-Moon rotating at $t_f$)', 'Interpreter', 'latex', 'FontSize', 8);
hold off;

% Results plot 1:
figure; hold on; grid on;

options_plot_LC = struct('color', 'r', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '-', 'DisplayName', 'LinCov 3$\sigma$ ellipse');
options_plot_UT = struct('color', 'b', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '--', 'DisplayName', 'UT 3$\sigma$ ellipse');
options_plot_MC = struct('color', 'g', 'alpha', 0.2, 'linewidth', 1.5, 'linestyle', '-', 'DisplayName', 'MC 3$\sigma$ ellipse');

plot(Sf_MC(1, :), Sf_MC(2, :), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.7, 0.7, 0.7], ...
     'MarkerEdgeColor', 'none', 'DisplayName', 'MC Samples at final time', 'MarkerSize', 2);          % Monte Carlo simulation points

plot_uncertainty_ellipse(mean_LC(1:2), cov_LC(1:2, 1:2), 3, options_plot_LC);                         % Ellipse of LinCov
plot_uncertainty_ellipse(mean_MC(1:2), cov_MC(1:2, 1:2), 3, options_plot_MC);                         % Ellipse of MC
plot_uncertainty_ellipse(mean_UT(1:2), cov_UT(1:2, 1:2), 3, options_plot_UT);                         % Ellipse of UT

plot(x_orbit, y_orbit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference final orbit');              % Target Orbit

plot(mean_MC(1), mean_MC(2), 'ro', 'MarkerSize', 14, 'LineWidth', 1.5, 'DisplayName', 'MC Mean');     % Mean of LinCov
plot(mean_MC(1), mean_MC(2), 'g+', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'MC Mean');     % Mean of MC
plot(mean_UT(1), mean_UT(2), 'bx', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'UT Mean');     % Mean of UT

xlabel('$x \, [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y \, [DU]$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
title('Zommed-in view of the comparison of Position Mean and Uncertainty Ellipses (@EMB Earth-Moon rotating at $t_f$)', 'Interpreter', 'latex', 'FontSize', 8);
axis([mean_UT(1)-0.00165 mean_UT(1)+0.00165 mean_UT(2)-0.00165 mean_UT(2)+0.00165])
hold off;

% Results plot 2
% Compute the recurring parameter of the time vector dimension:
m_t = length(t_vect);

% Preallocate matrices for eigenvalues
lambdas_posLC = zeros(2, m_t);
lambdas_posUT = zeros(2, m_t);
lambdas_posMC = zeros(2, m_t);

lambdas_velLC = zeros(2, m_t);
lambdas_velUT = zeros(2, m_t);
lambdas_velMC = zeros(2, m_t);

% Efficient computation
for j = 1:m_t
    lambdas_posLC(:, j) = eig(P_LC(1:2, 1:2, j));
    lambdas_posUT(:, j) = eig(P_UT(1:2, 1:2, j));
    lambdas_posMC(:, j) = eig(P_MC(1:2, 1:2, j));
    
    lambdas_velLC(:, j) = eig(P_LC(3:4, 3:4, j));
    lambdas_velUT(:, j) = eig(P_UT(3:4, 3:4, j));
    lambdas_velMC(:, j) = eig(P_MC(3:4, 3:4, j));
end

% Max eigenvalues and scaled results
te_posLC = 3 * sqrt(max(lambdas_posLC, [], 1));
te_posUT = 3 * sqrt(max(lambdas_posUT, [], 1));
te_posMC = 3 * sqrt(max(lambdas_posMC, [], 1));
te_velLC = 3 * sqrt(max(lambdas_velLC, [], 1));
te_velUT = 3 * sqrt(max(lambdas_velUT, [], 1));
te_velMC = 3 * sqrt(max(lambdas_velMC, [], 1));


% Define a consistent color scheme
colors = lines(3); % MATLAB's built-in line color palette

% Position evolution plot
figure
plot(t_vect, te_posLC, 'Color', colors(1, :), 'LineWidth', 1.5, ...
    'Marker','diamond', 'LineStyle', '-', 'DisplayName', 'LinCov $3\sqrt{\max{\lambda_i \left( \mathbf{P_{r}} \right)}}$');
hold on
plot(t_vect, te_posUT, 'Color', colors(2, :), 'LineWidth', 1.5, ...
    'LineStyle', '-', 'DisplayName', 'UT $3\sqrt{\max{\lambda_i \left( \mathbf{P_{r}} \right)}}$');
plot(t_vect, te_posMC, 'Color', colors(3, :), 'LineWidth', 1.5, ...
    'LineStyle', '--', 'DisplayName', 'MC $3\sqrt{\max{\lambda_i \left( \mathbf{P_{r}} \right)}}$');
hold off
xlabel('Time [TU]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$3\sqrt{\max{\lambda_i}}$ for Position', 'Interpreter', 'latex', 'FontSize', 14);
title('Time Evolution of Position Metrics (@EMB Earth-Moon rotating frame)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
grid on

% Velocity evolution plot
figure
plot(t_vect, te_velLC, 'Color', colors(1, :), 'LineWidth', 1.5, ...
    'Marker','diamond', 'LineStyle', '-', 'DisplayName', 'LinCov $3\sqrt{\max{\lambda_i \left( \mathbf{P_{v}} \right)}}$');
hold on
plot(t_vect, te_velUT, 'Color', colors(2, :), 'LineWidth', 1.5, ...
    'LineStyle', '-', 'DisplayName', 'UT $3\sqrt{\max{\lambda_i \left( \mathbf{P_{v}} \right)}}$');
plot(t_vect, te_velMC, 'Color', colors(3, :), 'LineWidth', 1.5, ...
    'LineStyle', '--', 'DisplayName', 'MC $3\sqrt{\max{\lambda_i \left( \mathbf{P_{v}} \right)}}$');
hold off
xlabel('Time [TU]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$3\sqrt{\max{\lambda_i}}$ for Velocity', 'Interpreter', 'latex', 'FontSize', 14);
title('Time Evolution of Velocity Metrics (@EMB Earth-Moon rotating frame)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
grid on


% Plot 3
n_components = size(Sf_MC, 1); % Extract the number of components in the state vector
figure; % Create a figure for the Q-Q plots
for i = 1:n_components
    % Subplot for each component
    subplot(ceil(n_components / 2), 2, i); % Arrange subplots in rows and columns

    % Q-Q plot for the i-th component of Monte Carlo samples
    h = qqplot(Sf_MC(i, :));
    
    % Set marker to points and reduce the marker size
    set(h(1), 'Marker', '+', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    
    % Set the title based on component number
    if i == 1
        comp_name = 'x-coordinate position';
    elseif i == 2
        comp_name = 'y-coordinate position';
    elseif i == 3
        comp_name = 'x-component velocity';
    elseif i == 4
        comp_name = 'y-component velocity';
    else
        comp_name = ['Component ' num2str(i)];
    end
    
    % Add title with LaTeX interpreter
    title(['Q-Q Plot for ', comp_name], 'Interpreter', 'latex', 'FontSize', 14);
    
    % Add labels with LaTeX interpreter
    xlabel('Theoretical Quantiles', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Sample Quantiles', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Enhance visual appeal: grid and axis limits
    grid on;
    axis tight;
    
    % Add legend with corrected items:
    legend([h(1), h(2)], {'Quantiles of the Sample Data', 'Theoretical Quantiles from Normal Dist.'}, ...
           'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
end

% Adjust layout for better visibility
sgtitle('Q-Q Plots for Monte Carlo Samples at Final Time (@EMB Earth-Moon rotating frame)', 'Interpreter', 'latex', 'FontSize', 10);

%% Functions

%--------------------------------------------------------------------------
function [dSdt] = xyPBRFBP_STM_ROT(t, S0, par)
% [dSdt] = xyPBRFBP_STM_ROT(t, S0, par)
% Description:
% Computes the right-hand side of the dynamics for the state vector under 
% the Planar Bicircular Restricted Four-Body Problem (PBRFBP) including 
% the evolution of the State Transition Matrix (STM).
%
% Inputs:
% t:   [1x1] Current time for the propagation [s]
% S0:  [20x1] State vector at the current time. The first four components 
%             represent the position and velocity in the planar PBRFBP 
%             system. The remaining 16 elements correspond to the 
%             reshaped 4x4 STM, stored as a column vector.
% par: [struct] Structure containing parameters for the system:
%       - mu:  [1x1] Gravitational parameter of the Earth-Moon system
%       - ms:  [1x1] Scaled mass of the Sun
%       - ws:  [1x1] Scaled angular velocity of the Sun
%       - rho: [1x1] Scaled distance between the Sun and Earth-Moon barycenter
%
% Outputs:
% dSdt: [20x1] Time derivative of the state vector. The first four 
%              components correspond to the derivatives of position and 
%              velocity. The remaining 16 elements correspond to the 
%              reshaped time derivative of the 4x4 STM, stored as a 
%              column vector.
     
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

    % Calculate the important radii:
    r1 = (mu + x)^2 + y^2;
    r2 = (mu + x - 1)^2 + y^2;
    r3 = (x - rho*cos(t*ws))^2 + (y - rho*sin(t*ws))^2;
    
    % Compute r1^(3/2), r2^(3/2), r3^(3/2) for more efficient calculations:
    r1_3_2 = r1^(3/2);
    r2_3_2 = r2^(3/2);
    r3_3_2 = r3^(3/2);  
    
    % Compute the derivative of the potential:
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
    dSdt(1:2) = S0(3:4);
    dSdt(3)   = dOM4dx + 2*vy;
    dSdt(4)   = dOM4dy - 2*vx;
    dSdt(5:end) = Phidot(:);
    
end

% -------------------------------------------------------------------------
function [Sf, PHIf, tf, xx, tt] = propagation(t0, S0, tf, par)

% [xf, PHIf, tf, xx, tt]  = propagation(t0, S0, tf, par)
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
% Outputs:
% Sf: [6x1] final state at tf
% PHIf: [6x6] STM at tf. This output is set to null matrix in case 
%       label='xyz_n_body' is chosen
% tf: [1x1] final time 
% xx: [6x?)] matrix containing the evolution of the state vector along the
%     integration
% tt: [1x??] vector containing the evolution of the time index along the
%     propagation
    
    % Initialize State Transition Matrix at t0
    Phi0 = eye(4);
    
    % Append to initial conditions the conditions for the STM
    x0Phi0 = [S0; Phi0(:)];
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', [1e-10*ones(1,2), 1e-13*ones(1,2), 1e-13*ones(1,16)]); % Set the options    
    [tt, xx] = ode78(@(t,x) xyPBRFBP_STM_ROT(t, x, par), [t0 tf], x0Phi0, options);                    % Integrate
    
    % Extract state vector and State Transition Matrix
    Sf = xx(end,1:4)';
    PHIf = reshape(xx(end, 5:end), 4, 4);
    tf = tt(end);

end

%--------------------------------------------------------------------------
function [E, P] = LinCov(S0, P0, t, par)
% [E, P] = LinCov(S0, P0, t, par)
% Description:
% Performs linear covariance propagation for uncertainty propagation using 
% the State Transition Matrix (STM). It computes the evolution of the mean 
% state and the covariance matrix over a given time vector.
%
% Inputs:
% S0:  [nx1] Initial state vector at the start of the propagation.
% P0:  [nxn] Initial covariance matrix at the start of the propagation.
% t:   [1xk] Time vector for propagation, where k is the number of time steps.
% par: [struct] Structure containing system parameters for propagation.
% propagation: [function handle] Nonlinear propagation function. It returns 
%              the propagated state and STM:
%              [Sf, PHIf, ~, ~, ~] = propagation(t_start, S0, t_end, par)
%
% Outputs:
% E: [nxk] Propagated mean state vector at each time step.
% P: [nxnxk] Propagated covariance matrix at each time step.
   
    % Dimension of the state vector
    n = length(S0);

    % Initialize Outputs
    E = zeros(n, length(t));    % Mean state at each time step
    P = zeros(n, n, length(t)); % Covariance matrix at each time step

    % Initialize mean and covariance
    E(:, 1) = S0;
    P(:, :, 1) = P0;

    % Loop through time steps
    for j = 2:length(t)

        % Propagate the solution (state and state transition matrix)
        [Sf, PHIf, ~, ~, ~] = propagation(t(1), S0, t(j), par);

        % Calculate the mean at the current time step
        E(:, j) = Sf;

        % Update the covariance matrix using the state transition matrix
        P(:, :, j) = PHIf*P0*PHIf';
    end
end

%--------------------------------------------------------------------------
function [E, P] = UnTr(S0, P0, t, alpha, beta, par)
% [E, P] = UnTr(S0, P0, t, alpha, beta, par)
% Performs uncertainty propagation using the Unscented Transform (UT).
% Inputs:
% - S0:    [nx1] Initial state mean vector.
% - P0:    [nxn] Initial covariance matrix.
% - t:     [1xm] Time vector.
% - alpha: [scalar] UT parameter controlling sigma point spread.
% - beta:  [scalar] UT parameter for prior distribution knowledge.
% - par:   [struct] Parameters for the propagation function.
% Outputs:
% - E: [nxm] Propagated mean state vector.
% - P: [nxnxm] Propagated covariance matrix.

    % Dimensions and scaling parameters
    n = length(S0);
    m_t = length(t);
    k = 0; % Secondary scaling parameter
    lambda = alpha^2 * (n + k) - n;
    n_lambda = n + lambda;

    % Weights for mean and covariance
    Wm = [lambda / n_lambda; repmat(1 / (2 * n_lambda), 2 * n, 1)];
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);

    % Initialization
    E = zeros(n, m_t);       % Mean state at each time step
    P = zeros(n, n, m_t);    % Covariance matrix at each time step
    E(:, 1) = S0;
    P(:, :, 1) = P0;

    % Scaled square root of initial covariance
    P_sqrt = sqrtm(n_lambda * P0);

    % Sigma points for t=1
    Chi = [S0, S0 + P_sqrt, S0 - P_sqrt];

    % Vectorize propagation for all time steps (except initial)
    for j = 2:m_t
        
        % Propagate all sigma points in parallel
        Chi_propagated = arrayfun(@(i) propagation(t(1), Chi(:, i), t(j), par), ...
                                  1:(2 * n + 1), 'UniformOutput', false);
        Chi_propagated = cat(2, Chi_propagated{:});

        % Compute propagated mean state
        E(:, j) = Chi_propagated * Wm;

        % Compute propagated covariance matrix
        diff = Chi_propagated - E(:, j);
        P(:, :, j) = diff * diag(Wc) * diff';
    end
end

%--------------------------------------------------------------------------
function plot_uncertainty_ellipse(mean, covariance, sig, options_plot)
% plot_uncertainty_ellipse(mean, covariance, sig, options_plot)
% Description:
% Plots a 2D uncertainty ellipse representing a Gaussian distribution
% with the given mean and covariance matrix, scaled by a specified 
% number of standard deviations.
%
% Inputs:
% mean:        [2x1] Vector representing the center of the ellipse.
% covariance:  [2x2] Covariance matrix defining the shape and orientation.
% sig:         [scalar] Scaling factor for the standard deviations (e.g., 1 for 1σ).
% options_plot: [struct] Structure specifying ellipse styling, containing:
%              - color: Line color (e.g., 'r', 'b', [0.5, 0.5, 0.5]).
%              - alpha: (Optional) Transparency level for the interior fill (0-1).
%              - linewidth: Line width of the ellipse edge.
%              - linestyle: Style of the edge line (e.g., '-', '--', ':').
%              - DisplayName: Name for the legend display.

    % Step 1: Eigenvalue decomposition for ellipse axes
    [V, D] = eig(covariance); % Eigenvectors (V) and eigenvalues (D)

    % Step 2: Sort eigenvalues (D) and corresponding eigenvectors (V)
    [D_sorted, idx] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order
    V = V(:, idx);                              % Reorder eigenvectors accordingly
    D_sorted = diag(D_sorted);                  % Reconstruct sorted eigenvalue matrix

    % Step 3: Generate unit circle points
    theta = linspace(0, 2*pi, 1000);            % Angle parameter for ellipse
    unit_circle = [cos(theta); sin(theta)];     % Points on a unit circle

    % Step 4: Scale and rotate the unit circle
    scaled_ellipse = sig * V * sqrt(D_sorted) * unit_circle;

    % Step 5: Translate the ellipse to the mean
    ellipse = bsxfun(@plus, scaled_ellipse, mean(:)); % Add mean to each ellipse point

    % Step 6: Plot the ellipse
    fill(ellipse(1, :), ellipse(2, :), options_plot.color, ...
        'EdgeColor', options_plot.color, ...         % Edge color of the ellipse
        'FaceAlpha', 0, ...                          % Fully transparent interior
        'LineWidth', options_plot.linewidth, ...     % Line width of the edge
        'LineStyle', options_plot.linestyle, ...     % Line style of the edge
        'DisplayName', options_plot.DisplayName);    % Name for legend display
end

%--------------------------------------------------------------------------
function [E, P, Sf_tf] = MC(S0, P0, t, par, n_samples)
% MC: Monte Carlo simulation for uncertainty propagation.
%
% Description:
% This function uses Monte Carlo sampling to propagate a state and
% covariance through a nonlinear system, returning the mean, covariance, 
% and final states at each time step.
%
% Inputs:
% S0:         [nx1] Initial state vector.
% P0:         [nxn] Initial covariance matrix.
% t:          [1xm] Time vector for propagation.
% par:        [struct] Parameters for the propagation function.
% n_samples:  [scalar] Number of Monte Carlo samples to use.
%
% Outputs:
% E:          [nxm] Mean state at each time step.
% P:          [nxnxm] Covariance matrix at each time step.
% Sf_tf:      [nxn_samples] Final sampled states at the last time step.

    % Step 1: Initialize useful quantities
    m_t = length(t);    % Number of time steps
    n = size(S0, 1);    % Dimension of the state vector
    t0 = t(1);          % Initial time for propagation

    % Step 2: Initialize output variables
    Sf = zeros(n, n_samples);  % Placeholder for final sampled states
    E = zeros(n, m_t);         % Placeholder for mean states
    P = zeros(n, n, m_t);      % Placeholder for covariance matrices

    % Step 3: Open parallel pool if not already open
    if isempty(gcp('nocreate'))
        parpool; % Start a default parallel pool
    end

    % Step 4: Monte Carlo simulation across time steps
    for j = 2:m_t
        tf = t(j); % Current final time to propagate to

        % Step 4.1: Perform parallel propagation of samples
        parfor i = 1:n_samples
            % Generate a random perturbation from the Gaussian distribution
            S_perturbed = mvnrnd(S0, P0);

            % Propagate the perturbed state (use row vector for `S_perturbed`)
            Sf(:, i) = propagation(t0, S_perturbed', tf, par);
        end

        % Step 4.2: Compute mean and covariance at current time step
        E(:, j) = mean(Sf, 2);     % Mean of sampled states
        P(:, :, j) = cov(Sf');     % Covariance matrix of sampled states

        % Step 4.3: Save sampled states at the final time step
        if j == m_t
            Sf_tf = Sf; % Store final sampled states
        end
    end
end
