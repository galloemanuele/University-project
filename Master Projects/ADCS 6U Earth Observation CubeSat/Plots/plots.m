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
%% PLOTS
% Run after the simulation from the main script

%% DISTURBING TORQUES
% Uncomment the Solar Radiation Pressure block from Environment for this
% plot
%
% M_GG = vecnorm(out.M_GG,2,2);
% M_SRP = vecnorm(out.M_SRP,2,2);
% M_B = vecnorm(out.M_B,2,2);
% 
% figure()
% semilogy(out.tout, M_GG, out.tout, M_SRP, out.tout, M_B,'Linewidth', 1.1)
% grid on
% yticks([1e-16 1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4]);
% xlabel('Time [s]','Interpreter','latex'); ylabel('Torque [Nm]','Interpreter','latex')
% legend('Gravity Gradient', 'Solar Radiation Pressure', 'Magnetic Field', 'Interpreter','latex', 'Location', 'southeast')

%% POINTING ERROR
figure()
plot(out.tout, out.pointing_err, 'LineWidth', 1)

xline(control.detumbling.MaxTime, '-.');
xline(control.slewManeuvre.MaxTime, '-.');
xline(out.tout(end), '-.');

txt1 = 'Detumbling';
t1 = text(100,70,txt1,'Interpreter','latex', 'FontSize', 10);
txt2 = 'Slew Manoeuvre';
t2 = text(1900,70,txt2,'Interpreter','latex','FontSize',10);
txt3 = 'Trajectory Tracking';
t3 = text(5500,70,txt3,'Interpreter','latex','FontSize',10);

xlabel('Time [s]','Interpreter','latex'); 
ylabel('Pointing Error [deg]', 'Interpreter','latex')

%% ANGULAR VELOCITY
figure()
plot(linspace(0, out.tout(end), length(out.omega_est_obs)), out.omega_est_obs,'LineWidth', 1) % With Control

xline(control.detumbling.MaxTime, '-.');
xline(control.slewManeuvre.MaxTime, '-.');
xline(out.tout(end), '-.');

xlabel('Time [s]','Interpreter','latex'); 
ylabel('Angular Velocity [rad/s]','Interpreter','latex');

txt1 = 'Detumbling';
t1 = text(100,3,txt1,'Interpreter','latex', 'FontSize', 14);
txt2 = 'Slew Manoeuvre';
t2 = text(1900,3,txt2,'Interpreter','latex','FontSize',14);
txt3 = 'Trajectory Tracking';
t3 = text(5500,3,txt3,'Interpreter','latex','FontSize',14);

xlabel('Time [s]','Interpreter','latex'); 
ylabel('Angular Velocity Uncontrolled [rad/s]','Interpreter','latex');

legend('$\omega_x$', '$\omega_y$','$\omega_z$', 'Interpreter','latex');

%% ACTUATORS TORQUE
figure()
plot(out.tout, out.MC2,'LineWidth', 1)

xline(control.detumbling.MaxTime, '-.');
xline(control.slewManeuvre.MaxTime, '-.');
xline(out.tout(end), '-.');

txt1 = 'Detumbling';
t1 = text(100,5e-3,txt1,'Interpreter','latex', 'FontSize', 14);
txt2 = 'Slew Manoeuvre';
t2 = text(1900,5e-3,txt2,'Interpreter','latex','FontSize',14);
txt3 = 'Trajectory Tracking';
t3 = text(5500,5e-3,txt3,'Interpreter','latex','FontSize',14);

xlabel('Time [s]','Interpreter','latex'); 
ylabel('Actuators Torque [Nm]','Interpreter','latex');
legend('$M_x$', '$M_y$','$M_z$', 'Interpreter','latex');

%% QUATERNION ERROR
figure()
plot(out.tout, (out.q_BN-out.q_LN))

xline(control.detumbling.MaxTime, '-.');
xline(control.slewManeuvre.MaxTime, '-.');
xline(out.tout(end), '-.');

txt1 = 'Detumbling';
t1 = text(270,0.1,txt1,'Interpreter','latex', 'FontSize', 14);
txt2 = 'Slew Manoeuvre';
t2 = text(2300,0.1,txt2,'Interpreter','latex','FontSize',14);
txt3 = 'Trajectory Tracking';
t3 = text(4600,0.1,txt3,'Interpreter','latex','FontSize',14);

xlabel('Time [s]','Interpreter','latex'); ylabel('Quaternion Error [-]', 'Interpreter','latex');
legend('$q_1$', '$q_2$','$q_3$', '$q_4$', 'Interpreter','latex', 'Location','southeast');

%% ANGULAR VELOCITY ERROR
figure()
plot(out.tout, (out.omega-out.w_LN_B))

xline(control.detumbling.MaxTime, '-.');
xline(control.slewManeuvre.MaxTime, '-.');
xline(out.tout(end), '-.');

txt1 = 'Detumbling';
t1 = text(270,0.075,txt1,'Interpreter','latex', 'FontSize', 14);
txt2 = 'Slew Manoeuvre';
t2 = text(2300,0.075,txt2,'Interpreter','latex','FontSize',14);
txt3 = 'Trajectory Tracking';
t3 = text(4600,0.075,txt3,'Interpreter','latex','FontSize',14);

xlabel('Time [s]','Interpreter','latex'); ylabel('Angular Velocity Error [rad/s]', 'Interpreter','latex');
legend('$\omega_1$', '$\omega_2$','$\omega_3$', 'Interpreter','latex', 'Location','southeast');

%% CONTROL FREE
% Control Logic and Actuators block should be commented for these plots
%% POINTING ERROR
% For this plot, the pointing error block has been transferred and is not
% present in the attached Simulink file

% figure()
% plot(out.tout, out.pointingerrornocontrol, 'LineWidth', 1) % Without control
% 
% xline(control.detumbling.MaxTime, '-.');
% xline(control.slewManeuvre.MaxTime, '-.');
% xline(out.tout(end), '-.');
% 
% txt1 = 'Detumbling';
% t1 = text(270,70,txt1,'Interpreter','latex', 'FontSize', 10);
% txt2 = 'Slew Manoeuvre';
% t2 = text(2300,70,txt2,'Interpreter','latex','FontSize',10);
% txt3 = 'Trajectory Tracking';
% t3 = text(4600,70,txt3,'Interpreter','latex','FontSize',10);
% 
% xlabel('Time [s]','Interpreter','latex'); 
% ylabel('Pointing Error [deg]', 'Interpreter','latex')

%% ANGULAR VELOCITY
% figure()
% plot(out.tout, out.omega, 'LineWidth', 1)
% 
% xline(control.detumbling.MaxTime, '-.');
% xline(control.slewManeuvre.MaxTime, '-.');
% xline(out.tout(end), '-.');
% 
% xlabel('Time [s]','Interpreter','latex'); 
% ylabel('Angular Velocity [rad/s]','Interpreter','latex');
% 
% txt1 = 'Detumbling';
% t1 = text(270,3,txt1,'Interpreter','latex', 'FontSize', 10);
% txt2 = 'Slew Manoeuvre';
% t2 = text(2300,3,txt2,'Interpreter','latex','FontSize',10);
% txt3 = 'Trajectory Tracking';
% t3 = text(4600,3,txt3,'Interpreter','latex','FontSize',10);
% 
% xlabel('Time [s]','Interpreter','latex'); 
% ylabel('Angular Velocity Uncontrolled [rad/s]','Interpreter','latex');
% 
% legend('$\omega_x$', '$\omega_y$','$\omega_z$', 'Interpreter','latex');

%% RANDOM ANALYSIS
%% Simulation
% Tsample = 0.1;
% for i = 1:5
%     w0 = 5*rand(3,1);
%     out(i) = sim('Project26.slx', 'FixedStep', '0.1', 'Solver', 'ode4', 'StopTime', '2*pi/n');
% end

%% POINTING ERROR RANDOM ANALYSIS
% close all
% figure(1)
% for i = 1:5
%     plot(out(i).tout, out(i).pointing_err, 'LineWidth', 1)
%     hold on
% 
% end
% 
% legend('p(w_1)', 'p(w_2)', 'p(w_3)', 'p(w_4)', 'p(w_5)')
% hold off
% xline(control.detumbling.MaxTime, '-.');
% xline(control.slewManeuvre.MaxTime, '-.');
% xline(out.tout(end), '-.');
% 
% txt1 = 'Detumbling';
% t1 = text(270,70,txt1,'Interpreter','latex', 'FontSize', 14);
% txt2 = 'Slew Manoeuvre';
% t2 = text(2300,70,txt2,'Interpreter','latex','FontSize',14);
% txt3 = 'Trajectory Tracking';
% t3 = text(4600,70,txt3,'Interpreter','latex','FontSize',14);
% 
% xlabel('Time [s]','Interpreter','latex'); 
% ylabel('Pointing Error - Random Analysis [deg]', 'Interpreter','latex')

%% ANGULAR VELOCITY RANDOM ANALYSIS
% The norm of w observed has been separately computed on Simulink and it's
% not present in the attached file

% close all
% figure(2)
% for i = 1:5
%     plot(linspace(0, out(i).tout(end), length(out(i).wobs_norm)), out(i).wobs_norm,'LineWidth', 1)
%     hold on
% end
% 
% legend('$\omega_1$', '$\omega_2$','$\omega_3$','$\omega_4$','$\omega_5$', 'Interpreter','latex');
% xline(control.detumbling.MaxTime, '-.');
% xline(control.slewManeuvre.MaxTime, '-.');
% xline(out.tout(end), '-.');
% 
% xlabel('Time [s]','Interpreter','latex'); 
% ylabel('Angular Velocity [rad/s]','Interpreter','latex');
% 
% txt1 = 'Detumbling';
% t1 = text(270,3,txt1,'Interpreter','latex', 'FontSize', 14);
% txt2 = 'Slew Manoeuvre';
% t2 = text(2300,3,txt2,'Interpreter','latex','FontSize',14);
% txt3 = 'Trajectory Tracking';
% t3 = text(4600,3,txt3,'Interpreter','latex','FontSize',14);
