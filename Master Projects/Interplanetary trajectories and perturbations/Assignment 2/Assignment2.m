%% Main script for Planetary Explorer Mission (Assignment 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       ORBITAL MECHANICS                                 %
%                    Academic year 2023/2024                              %
%                    M.Sc. Space Engineering                              %
%                     Politecnico di Milano                               %
%                                                                         %
%              Planetary Explorer Mission assignment                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Group N.2337
%
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta

clear 
close all
clc

% Path to the functions:
addpath(strcat(pwd,'/functions'));
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% Define the nominal orbit:
a_n=4.2164e+4;                          % semi-major axis [km]
e_n=0.0007;                             % eccentricity [-]
i_n=deg2rad(0.1378);                    % inclination [rad]
OM_n=deg2rad(0);                        % right ascension of the ascending node [rad]
w_n=deg2rad(0);                         % argument of perigee [rad]
f0_n=deg2rad(40);                       % initial true anomaly [rad]
par_n=[a_n, e_n, i_n, OM_n, w_n, f0_n]; % nominal keplerian elements

% Other useful parameters:
mu_E = astroConstants(13);  % gravitational parameter of the Earth [km^3/s^2]
mu_M = astroConstants(20);  % gravitational parameter of the Moon [km^3/s^2]
J2 = astroConstants(9);     % second zonal harmonic of the Earth
R_E = astroConstants(23);   % mean equatorial radius of the Earth [km]
T=2*pi*sqrt(a_n^3/mu_E);    % orbital period

%% Task selection:
% Running this section, you can select the task you want to perform:
Perform_task='true';

while strcmp(Perform_task, 'true')==1 

    task = input(['Choose the value of the variable "task" depending on what you want to be \n ' ...
                  'computed and then you can run only this section:\n' ...
                  '1 : Nominal and Repeating Ground Track (Perturbed and Unperturbed case)\n' ...
                  '2 : Orbit propagation with perturbations (Gauss Planetary equations \n    and Cartesian methods accounting for J2 and moon perturbations)\n    and comparison between the algoritmhs\n' ...
                  '3 : (Animated) Orbit evolution 3D representation \n' ...
                  '4 : (Non animated) Orbit evolution 3D representation \n' ...
                  '5 : Filtering of keplerian elements \n' ...
                  '6 : Comparison with real satellite data \n'...
                  '7 : None\n'...
                  'The task chosen is: \n']);
    
    switch task
        case 1
            %% Nominal Ground Track (Perturbed and Unperturbed case)
            
            % Determination of meaningful time intervals:
            N=10000;                                         % number of points for the timespan discretization
            t0_GT = 0;                                       % initial time [s]
            t_span1=linspace(t0_GT,T,N);                     % timespan of 1 orbit [s]
            t_span2=linspace(t0_GT,3600*24*365.25,N);        % timespan of 1 week [s]
            t_span3=linspace(t0_GT,3600*24*365.25*5,N);      % timespan of 30 days [s]
            
            
            thetaG0=0;                  % right ascension of Greenwich meridian at t0 [rad]
            w_E=deg2rad(15.04/3600);    % earth angular rate [rad/s]
            
            % Timespan 1: 
            [~,~,lon1,lat1]=groundTrack("kep",par_n,thetaG0,t_span1,mu_E);                                   % unperturbed 
            [~,~,lon1p,lat1p]=groundTrack_J2_moon("kep", par_n, thetaG0, t_span1, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,1)
            plot(lon1, lat1,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon1p, lat1p,'--k','LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon1(1), lat1(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon1(1),lat1(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon1(end), lat1(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon1(end),lat1(end), "End (Unperturbed)", "Color",'c')
            plot(lon1p(end), lat1p(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon1p(end),lat1p(end), "End (Perturbed)", "Color",'c')
            title('\bf{One orbit}', fontsize=12)
            %legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case')
            hold off
            
            % Timespan 2: 
            [~,~,lon2,lat2]=groundTrack("kep",par_n,thetaG0,t_span2,mu_E);                                   % unperturbed 
            [~,~,lon2p,lat2p]=groundTrack_J2_moon("kep", par_n, thetaG0, t_span2, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,2)
            plot(lon2, lat2,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon2p, lat2p,'--k', 'LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon2(1), lat2(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon2(1),lat2(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon2(end), lat2(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon2(end),lat2(end), "End (Unperturbed)", "Color",'c')
            plot(lon2p(end), lat2p(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon2p(end),lat2p(end), "End (Perturbed)", "Color",'c')
            title('\bf{One week}', fontsize=12)
            %legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case')
            hold off
            
            % Timespan 3: 
            [~,~,lon3,lat3]=groundTrack("kep",par_n,thetaG0,t_span3,mu_E);                                   % unperturbed 
            [~,~,lon3p,lat3p]=groundTrack_J2_moon("kep", par_n, thetaG0, t_span3, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,3)
            plot(lon3, lat3,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon3p, lat3p,'--k', 'LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon3(1), lat3(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon3(1),lat3(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon3(end), lat3(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon3(end),lat3(end), "End (Unperturbed)", "Color",'c')
            plot(lon3p(end), lat3p(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon3p(end),lat3p(end), "End (Perturbed)", "Color",'c')
            title('\bf{30 days}', fontsize=12)
            legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case', fontsize=10)
            hold off
            
            %% Repeating Ground Track (Unperturbed and Perturbed case):
            k=1;  % number of revolutions of the s/c to obtain the ground track repetition
            m=1;  % number of revolutions of the planet to obtain the ground track repetition
            
            % Modify the semimajor axis to obtain a repeating ground track:
            a_rep=a_GroundTrack_rep('unperturbed', m, k, mu_E, e_n, i_n, J2, R_E, w_E); % semimajor axis to obtain a repeating ground track in the unperturbed case
            a_rep_p=a_GroundTrack_rep('perturbed', m, k, mu_E, e_n, i_n, J2, R_E, w_E); % semimajor axis to obtain a repeating ground track in the perturbed case
            par_rep=[a_rep, par_n(2:6)];                                                % modified keplerian elements for unperturbed case
            par_rep_p=[a_rep_p, par_n(2:6)];                                            % modified keplerian elements for perturbed case
            T_rep=2*pi*sqrt((a_rep^3)/mu_E);                                            % orbit period related to the semimajor axis to obtain a repeating ground track
            T_rep_p=2*pi*sqrt((a_rep_p^3)/mu_E);                                        % orbit period related to the semimajor axis to obtain a repeating ground track
            
            % Definition of relevant time span:
            t_span_rep1=linspace(t0_GT, T_rep, N);          % one orbit for with the repeated ground track
            t_span_rep1_p=linspace(t0_GT, T_rep_p, N);      % one orbit for with the repeated ground track for the perturbed case
            t_span_rep2=linspace(t0_GT,3600*24*365.25,N);        % timespan of 1 week [s]
            t_span_rep3=linspace(t0_GT,3600*24*365.25*5,N);       % timespan of 1 month [s]
            
            % Timespan 1: 
            [~,~,lon1_rep,lat1_rep]=groundTrack("kep",par_rep,thetaG0,t_span_rep1,mu_E);                                       % unperturbed 
            [~,~,lon1p_rep,lat1p_rep]=groundTrack_J2_moon("kep", par_rep_p, thetaG0, t_span_rep1_p, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,4)
            plot(lon1_rep, lat1_rep,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon1p_rep, lat1p_rep,'--k', 'LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon1_rep(1), lat1_rep(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon1_rep(1),lat1_rep(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon1_rep(end), lat1_rep(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon1_rep(end),lat1_rep(end), "End (Unperturbed)", "Color",'c')
            plot(lon1p_rep(end), lat1p_rep(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon1p_rep(end),lat1p_rep(end), "End (Perturbed)", "Color",'c')
            title('\bf{Modified for one orbit}', fontsize=12)
            %legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case')
            hold off
            
            % Timespan 2: 
            [~,~,lon2_rep,lat2_rep]=groundTrack("kep",par_rep,thetaG0,t_span_rep2,mu_E);                                     % unperturbed 
            [~,~,lon2p_rep,lat2p_rep]=groundTrack_J2_moon("kep", par_rep_p, thetaG0, t_span_rep2, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,5)
            plot(lon2_rep, lat2_rep,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon2p_rep, lat2p_rep,'--k', 'LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon2_rep(1), lat2_rep(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon2_rep(1),lat2_rep(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon2_rep(end), lat2_rep(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon2_rep(end),lat2_rep(end), "End (Unperturbed)", "Color",'c')
            plot(lon2p_rep(end), lat2p_rep(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon2p_rep(end),lat2p_rep(end), "End (Perturbed)", "Color",'c')
            title('\bf{Modified for one week}', fontsize=12)
            %legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case')
            hold off
            
            % Timespan 3: 
            [~,~,lon3_rep,lat3_rep]=groundTrack("kep",par_rep,thetaG0,t_span_rep3,mu_E);                                     % unperturbed 
            [~,~,lon3p_rep,lat3p_rep]=groundTrack_J2_moon("kep", par_rep_p, thetaG0, t_span_rep3, mu_E, J2, mu_M, R_E, w_E); % perturbed
            subplot(2,3,6)
            plot(lon3_rep, lat3_rep,'LineStyle','none','Marker','.', 'LineWidth',0.5, 'Color','r')
            hold on
            plot(lon3p_rep, lat3p_rep,'--k', 'LineWidth',0.5)
            xlim([-180,180])
            ylim([-90,90])
            I=imread('EarthTexture.jpg');
            X=[-180,180];
            Y=[90,-90];
            h=image(X,Y,I);
            uistack(h,'bottom')
            plot(lon3_rep(1), lat3_rep(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon3_rep(1),lat3_rep(1), "Start (Unperturbed and Perturbed)", "Color",'c')
            plot(lon3_rep(end), lat3_rep(end), 'LineWidth', 2,'Color','c', 'Marker','o' )
            text(lon3_rep(end),lat3_rep(end), "End (Unperturbed)", "Color",'c')
            plot(lon3p_rep(1), lat3p_rep(1), 'LineWidth', 2,'Color','c', 'Marker','o')
            text(lon3p_rep(end),lat3p_rep(end), "End (Perturbed)", "Color",'c')
            title('\bf{Modified for 30 days}', fontsize=12)
            %legend('Unperturbed case', 'Perturbed case', 'Initial point of unperturbed and perturbed case', 'Final point of unperturbed case', 'Final point of perturbed case')
            hold off
    
            % Algorithm to choose another task:
            choice1=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice1==2
                Perform_task='false';
            end 
        case 2
            %% Propagate the orbit with the J2 and moon perturbations
            % Choose an adeguate time span for the propagations:
            date0=[2028, 01, 01, 00, 00, 00];           % initial date for propagation
            t0=date2mjd2000(date0)*60*60*24;            % initial date for propagation [s]
            k_prop=500;                                 % number of orbits chosen for the propagation
            N_prop=10000;                               % chosen discretization for the propagation
            tspan=linspace(t0, t0+k_prop*T ,N_prop);    % timespan of the discretization
            
            % Set options for the ODE solver:
            options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
            
            % Propagation through Cartesian coordinates:
            [r0,v0]=kep2car(a_n, e_n, i_n, OM_n, w_n, f0_n, mu_E);
            y0=[r0; v0]; % initial condition for the propagation
            tic          % useful to calculate the computational time used by the algorithm
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'XYZ');
            dy_car=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'XYZ');
            [T_Car, Y_car]=ode113(dy_car, tspan, y0, options);
            Computational_Time_car=toc;
            acar = zeros(length(tspan),1);
            ecar = zeros(length(tspan),1);
            icar = zeros(length(tspan),1);
            OMcar = zeros(length(tspan),1);
            wcar = zeros(length(tspan),1);
            fcar = zeros(length(tspan),1);
            for i = 1:size(Y_car,1)
                r = Y_car(i, 1:3);
                v = Y_car(i, 4:6);
                [acar(i), ecar(i), icar(i), OMcar(i), wcar(i), fcar(i)] = car2kep(r', v', mu_E);
            end
            OMcar = unwrap(OMcar);
            wcar = unwrap(wcar-2*pi);
            fcar = unwrap(fcar);
            
            % Propagation through Gauss's planetary equations, whose initial condition are the ones of nominal orbit:
            tic
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'RSW');
            dy_Gauss=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'RSW');
            [T_Gauss, Y_Gauss]=ode113(dy_Gauss, tspan, par_n, options);
            Computational_Time_Gauss=toc;
            aGauss = Y_Gauss(:,1);
            eGauss = Y_Gauss(:,2);
            iGauss = Y_Gauss(:,3);
            OMGauss = unwrap(Y_Gauss(:,4));
            wGauss = unwrap(Y_Gauss(:,5));
            fGauss = unwrap(Y_Gauss(:,6));
            
            
            %% Plot the history of the Keplerian elements:
            % Errors between Gauss and Cartesian propagation:
            err_a_plot = abs(acar-aGauss)./abs(a_n);          % relative error of semi-major axis
            err_e_plot = abs(ecar-eGauss);                    % absolute error of eccentricity
            err_i_plot = abs(icar-iGauss)/(2*pi);             % relative error of inclination
            err_OM_plot = abs(OMcar-OMGauss)/(2*pi);          % relative error of right ascension of ascending node
            err_w_plot = abs(wcar-wGauss)/(2*pi);             % relative error of argoument of periapsis
            err_th_plot = abs(fcar-fGauss)./abs(fGauss);      % relative error of true anomaly
            
            % Plot of comparison between keplerian elements with Gauss and
            % Cartesian propagation:
            figure()
            subplot(2,3,1)
            plot((T_Car-t0)/T,acar,'LineWidth',0.5)
            title('\bf{Evolution of $a$ Cartesian method [km]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$a_{car}$ [km]', fontsize=12);
            grid on
            
            subplot(2,3,2)
            plot((T_Gauss-t0)/T, aGauss,'LineWidth',0.5)
            title('\bf{Evolution of $a$  Gauss method [km]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$a_{Gauss}$ [km]', fontsize=12);
            grid on
            
            subplot(2,3,3)
            semilogy((T_Gauss-t0)/T,err_a_plot,'LineWidth',0.5)
            title('\bf{Relative error of $a$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|a_{car} - a_{Gauss}|/|a_0|$', fontsize=12);
            grid on
            
            subplot(2,3,4)
            plot((T_Car-t0)/T,ecar,'LineWidth',0.5)
            title('\bf{Evolution of $e$ Cartesian method [-]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$e_{car} [-]$', fontsize=12);
            grid on
            
            subplot(2,3,5)
            plot((T_Gauss-t0)/T,eGauss,'LineWidth',0.5)
            title('\bf{Evolution of $e$ Gauss method [-]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$e_{Gauss} [-]$', fontsize=12);
            grid on
            
            subplot(2,3,6)
            semilogy((T_Gauss-t0)/T,err_e_plot,'LineWidth',0.5)
            title('\bf{Absolute error of $e$ [-]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|e_{Car} - e_{Gauss}| [-]$', fontsize=12);
            grid on
            
            figure()
            subplot(2,3,1)
            plot((T_Car-t0)/T,rad2deg(icar),'LineWidth',0.5)
            title('\bf{Evolution of $i$ Cartesian method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$i_{car} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,2)
            plot((T_Gauss-t0)/T,rad2deg(iGauss),'LineWidth',0.5)
            title('\bf{Evolution of $i$ Gauss method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$i_{Gauss} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,3)
            semilogy((T_Gauss-t0)/T,rad2deg(err_i_plot),'LineWidth',0.5)
            title('\bf{Relative error of $i$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|i_{car} - i_{gauss}|/360^{\circ}$', fontsize=12);
            grid on
            
            subplot(2,3,4)
            plot((T_Car-t0)/T,rad2deg(OMcar),'LineWidth',0.5)
            title('\bf{Evolution of $\Omega$ Cartesian method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\Omega_{car} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,5)
            plot((T_Gauss-t0)/T,rad2deg(OMGauss),'LineWidth',0.5)
            title('\bf{Evolution of $\Omega$ Gauss method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\Omega_{Gauss} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,6)
            semilogy((T_Gauss-t0)/T,rad2deg(err_OM_plot),'LineWidth',0.5)
            title('\bf{Relative error of $\Omega$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\Omega_{car} - \Omega_{Gauss}|/360^{\circ}$', fontsize=12);
            grid on
            
            figure()
            subplot(2,3,1)
            plot((T_Car-t0)/T,rad2deg(wcar),'LineWidth',0.5)
            title('\bf{Evolution of $\omega$ Cartesian method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\omega_{car} [deg]$', fontsize=12);
            grid on
                        
            subplot(2,3,2)
            plot((T_Gauss-t0)/T,rad2deg(wGauss),'LineWidth',0.5)
            title('\bf{Evolution of $\omega$ Gauss method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\omega_{Gauss} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,3)
            semilogy((T_Gauss-t0)/T,rad2deg(err_w_plot),'LineWidth',0.5)
            title('\bf{Relative error of $\omega$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\omega_{car} - \omega_{Gauss}|/360^{\circ}$', fontsize=12);
            grid on
            
            subplot(2,3,4)
            plot((T_Car-t0)/T,rad2deg(fcar),'LineWidth',0.5)
            title('\bf{Evolution of $\theta$ Cartesian method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\theta_{car} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,5)
            plot((T_Gauss-t0)/T,rad2deg(fGauss),'LineWidth',0.5)
            title('\bf{Evolution of $\theta$ Gauss method [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\theta_{Gauss} [deg]$', fontsize=12);
            grid on
            
            subplot(2,3,6)
            semilogy((T_Gauss-t0)/T,rad2deg(err_th_plot),'LineWidth',0.5)
            title('\bf{Relative error of $\theta$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\theta_{car} - \theta_{Gauss}|/360^{\circ}$', fontsize=12);
            grid on
            hold off
            
            %% Comparison between the two algoritmhs in terms of time 
            
            % Comparison for the specific case analysed:
            if Computational_Time_car>Computational_Time_Gauss
                TimeSave=Computational_Time_car-Computational_Time_Gauss;
                fprintf('Propagation by means of Gauss''s planetary equations \n for the case analysed is more \n time efficient of TimeSave=%f s .\n',TimeSave)
            elseif Computational_Time_car<Computational_Time_Gauss
                TimeSave=Computational_Time_Gauss-Computational_Time_car;
                fprintf('Propagation by means of Cartesian coordinates \n for the case analysed is more \n time efficient of TimeSave=%f s .\n',TimeSave)
            else
                fprintf('Propagation by means of Gauss''s planetary equations \n and Cartesian coordinates are equivalent for the case analysed \n with respect of time consumption .\n')
            end

            % General comparison of the two algoritmhs in terms of time:
            t0_comp=0;                                                      % initial time for plotting orbit
            k_comp=100;                                                     % number of orbit propagation for the evolution plot
            tspan_max=k_comp*T;                                             % Maximum time span for the integration
            n_points=k_comp*linspace(50, 1000, 50);                         % Array of number of time steps for each integration
            
            % Preallocate space for the vectors of time consumption for the two methods:
            t_car_comp=zeros(1,length(n_points));
            t_Gauss_comp=zeros(1,length(n_points));

            fprintf('Calculating the comparison between Gauss''s planetary equations and Cartesian method... \n ')

            % Cycle over the different time spans to compare the methods in terms of time:
            for i=1:length(n_points)
                
                % For each loop iteration change the points considered in the time span:
                tspan_comp=linspace(t0_comp, tspan_max, n_points(i));

                % Cartesian method (initial condition already set for the case analysed):
                tic                                                        % useful to calculate the computational time used by the algorithm
                acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'XYZ');
                dy_car=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'XYZ');
                [T_Car_comp, Y_car_comp] = ode113(dy_car, tspan_comp, y0, options);
                t_car_comp(i) = toc;
                
                % Gauss's planetary equations (initial conditions already set):
                tic
                acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'RSW');
                dy_Gauss=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'RSW');
                [T_Gauss_comp, Y_Gauss_comp]=ode113(dy_Gauss, tspan_comp, par_n, options);
                t_Gauss_comp(i) = toc;

            end

            fprintf('Finished...')

            % Plot:
            figure()
            semilogy(n_points, t_car_comp)
            hold on
            semilogy(n_points, t_Gauss_comp)
            legend('Cartesian method','Gauss''s planetary equations method'),
            title('Computational time'),
            xlabel('Number of time steps'),
            ylabel('Time [s]'),
            xlim([n_points(1), n_points(end)])
            ylim([0, 2.5])
            grid on


            % Algorithm to choose another task:
            choice2=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice2==2
                Perform_task='false';
            end 

        case 3
            %% (Animated) Perturbed Orbit Evolution 3D Representation
            % Calculating perturbed parameters of orbits:
            tic
            fprintf('Calculating the perturbed orbits...\n')
            t0_3D_plot=0;                                                      % initial time for plotting orbit
            k_evo=50000;                                                       % number of orbit propagation for the evolution plot
            N_evo=40;                                                          % number of discretization for the propagation for the evolution plot for a single orbit
            tspan_evo = linspace(t0_3D_plot,t0_3D_plot+k_evo*T ,N_evo*k_evo);  % timespan for propagation for the evolution [s]
    
            [r0_plot, v0_plot]=kep2car(par_n(1), par_n(2), par_n(3),par_n(4), par_n(5), par_n(6), mu_E);
            y0_plot=[r0_plot; v0_plot];
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'XYZ');
            dy_Gauss=@(t, y) eq_motion(t, y, acc_pert(t,y), mu_E, 'XYZ');
            [T_Gauss_car_plot3D, Y_car_plot3D]=ode113(dy_Gauss, tspan_evo, y0_plot, options);  % propagation of the orbit through Gauss equations
            toc
            fprintf('Finished\n')
    
            %% Plot:
            fprintf('Orbit evolution representation calculations in progress...\n');
            tic
            
            fig=figure();
            % Maximize the figure to full screen
            set(fig, 'Position', get(0, 'Screensize'));
    
            % Create an annotation for the progress bar
            progress_text = annotation('textbox', [0.01, 0.9, 0.1, 0.1], ...
                'String', 'Orbit plotting: 0', ...
                'FitBoxToText', 'on', ...
                'EdgeColor', 'none', ...
                'BackgroundColor', 'white', ...
                'FontSize', 14, ...
                'FontWeight', 'bold', ...
                'Color', 'blue');
    
            % Set color axis limits for the colorbar
            clim([1, k_evo]);
    
            % Create a custom colormap from red to blue
            colormap_red_to_blue = [linspace(1, 0, k_evo); zeros(1, k_evo); linspace(0, 1, k_evo)]';
                 
            % Initialize color data for the colorbar
            color_data = zeros(size(Y_car_plot3D, 1), 1);
           
            % Plot the Earth
            Earth3D
            hold on
            title('Perturbed orbits for J2 and Moon (animated)');
            grid on
            xlabel('x [km]')
            ylabel('y [km]')
            zlabel('z [km]')
            xlim([-50000, 50000])
            ylim([-50000, 50000])
            zlim([-50000, 50000])
    
            % Set the desired view angles:
            azimuth = 45;             % azimuth angle (degrees)
            elevation = 30;           % elevation angle (degrees)        
            view(azimuth, elevation); % set the 3D view
   
    
            for i=[1:1000:k_evo, k_evo] % Plot one orbit per cycle
                if i==1 % select the first orbit
                   % Plot the first orbit:
                   plot3(Y_car_plot3D(1:N_evo*i, 1), Y_car_plot3D(1:N_evo*i, 2), Y_car_plot3D(1:N_evo*i, 3), '--r','LineWidth', 3);
                   drawnow;
                elseif i==k_evo % select the last orbit
                    % Plot the last orbit
                    plot3(Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 1), Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 2), Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 3), '--b','LineWidth', 3);
                    drawnow;
                else % select the other orbits
                    % Plot the intermediate orbits:
                    pib=plot3(Y_car_plot3D(N_evo*(i-1):N_evo*i, 1), Y_car_plot3D(N_evo*(i-1):N_evo*i, 2), Y_car_plot3D(N_evo*(i-1):N_evo*i, 3),'LineWidth', 0.5, 'Color', colormap_red_to_blue(i,:));
                    pib.HandleVisibility = 'off'; % delete these plot from the legend visibility
                    drawnow;
                end 
    
                % Add a pause to control the animation speed (adjust the duration as needed)
                pause(0.05);
    
                % Update color data for the colorbar
                color_data(i) = min(i, k_evo);
      
                % Update progress bar text
                set(progress_text, 'String', sprintf('Orbit plotting: %d', i));
            end 
                
            % Plot the initial point:
            plot3(Y_car_plot3D(1, 1), Y_car_plot3D(1, 2), Y_car_plot3D(1, 3), 'Marker','o','color', 'red', 'LineWidth', 3);
            % Plot the final point:
            plot3(Y_car_plot3D(end, 1), Y_car_plot3D(end, 2), Y_car_plot3D(end, 3), 'Marker','o','color', 'blue', 'LineWidth', 3);
    
            % Build the legend
            legend('', 'Initial Orbit','Final Orbit','Starting Point', 'Ending Point')
            toc
           
            % Create custom colorbar with specific tick values based on color_data
            colorbar_handle = colorbar;
            desired_ticks = 0:5000:k_evo;
            set(colorbar_handle, 'Ticks', linspace(-7473, 5731, length(desired_ticks)), 'TickLabels', desired_ticks, "Limits", [-7473, 5731], 'ColorMap', colormap_red_to_blue(1:end,:));
            colorbar_title = get(colorbar_handle, 'Title');
            set(colorbar_title,'interpreter','latex', 'String', 'Orbit Number');
            
            fprintf('Finished.\n');
            
            % Algorithm to choose another task:
            choice3=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice3==2
                Perform_task='false';
            end 
        case 4
            %% (Non animated) Pertubed Orbit Evolution 3D Representation
            % Calculating perturbed parameters of orbits:
            tic
            fprintf('Calculating the perturbed orbits...\n')
            t0_3D_plot=0;                                                      % initial time for plotting orbit
            k_evo=50000;                                                       % number of orbit propagation for the evolution plot
            N_evo=40;                                                          % number of discretization for the propagation for the evolution plot for a single orbit
            tspan_evo = linspace(t0_3D_plot,t0_3D_plot+k_evo*T ,N_evo*k_evo);  % timespan for propagation for the evolution [s]
    
            [r0_plot, v0_plot]=kep2car(par_n(1), par_n(2), par_n(3),par_n(4), par_n(5), par_n(6), mu_E);
            y0_plot=[r0_plot; v0_plot];
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'XYZ');
            dy_Gauss=@(t, y) eq_motion(t, y, acc_pert(t,y), mu_E, 'XYZ');
            [T_Gauss_car_plot3D, Y_car_plot3D]=ode113(dy_Gauss, tspan_evo, y0_plot, options);  % propagation of the orbit through Gauss equations
            toc
            fprintf('Finished\n')
    
            % Plot:
            fprintf('Orbit evolution representation calculations in progress...\n');
            tic
            
            fig=figure();
            % Maximize the figure to full screen
            set(fig, 'Position', get(0, 'Screensize'));
            % Set color axis limits for the colorbar
            clim([1, k_evo]);
    
            % Create a custom colormap from red to blue
            colormap_red_to_blue = [linspace(1, 0, k_evo); zeros(1, k_evo); linspace(0, 1, k_evo)]';
                 
            % Initialize color data for the colorbar
            color_data = zeros(size(Y_car_plot3D, 1), 1);
           
            % Plot the Earth
            Earth3D
            hold on
            title('\bf{Perturbed orbits for J2 and Moon (animated)}', fontsize=12);
            grid on
            xlabel('x [km]', fontsize=12)
            ylabel('y [km]', fontsize=12)
            zlabel('z [km]', fontsize=12)
            xlim([-50000, 50000])
            ylim([-50000, 50000])
            zlim([-50000, 50000])
    
            % Set the desired view angles:
            azimuth = 45;             % azimuth angle (degrees)
            elevation = 30;           % elevation angle (degrees)        
            view(azimuth, elevation); % set the 3D view
            for i=[1:1000:k_evo, k_evo] % Plot one orbit per cycle
                if i==1 % select the first orbit
                   % Plot the first orbit:
                   plot3(Y_car_plot3D(1:N_evo*i, 1), Y_car_plot3D(1:N_evo*i, 2), Y_car_plot3D(1:N_evo*i, 3), '--r','LineWidth', 3);
                elseif i==k_evo % select the last orbit
                    % Plot the last orbit
                    plot3(Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 1), Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 2), Y_car_plot3D(N_evo*(k_evo-1):N_evo*k_evo, 3), '--b','LineWidth', 3);
                else % select the other orbits
                    % Plot the intermediate orbits:
                    pib=plot3(Y_car_plot3D(N_evo*(i-1):N_evo*i, 1), Y_car_plot3D(N_evo*(i-1):N_evo*i, 2), Y_car_plot3D(N_evo*(i-1):N_evo*i, 3),'LineWidth', 0.5, 'Color', colormap_red_to_blue(i,:));
                    pib.HandleVisibility = 'off'; % delete these plot from the legend visibility
                end 
            
                % Update color data for the colorbar
                color_data(i) = min(i, k_evo);
      
            end 
                
            % Plot the initial point:
            plot3(Y_car_plot3D(1, 1), Y_car_plot3D(1, 2), Y_car_plot3D(1, 3), 'Marker','o','color', 'red', 'LineWidth', 3);
            % Plot the final point:
            plot3(Y_car_plot3D(end, 1), Y_car_plot3D(end, 2), Y_car_plot3D(end, 3), 'Marker','o','color', 'blue', 'LineWidth', 3);
    
            % Build the legend
            legend('', 'Initial Orbit','Final Orbit','Starting Point', 'Ending Point')
            toc
    
            % Create custom colorbar with specific tick values based on color_data
            colorbar_handle = colorbar;
            desired_ticks = 0:5000:k_evo;
            set(colorbar_handle, 'Ticks', linspace(-7473, 5731, length(desired_ticks)), 'TickLabels', desired_ticks, "Limits", [-7473, 5731], 'ColorMap', colormap_red_to_blue(1:end,:));
            colorbar_title = get(colorbar_handle, 'Title');
            set(colorbar_title,'interpreter','latex', 'String', 'Orbit Number');
            
            fprintf('Finished.\n');

            % Algorithm to choose another task:
            choice4=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice4==2
                Perform_task='false';
            end 
        case 5
            %% Filtering of high frequencies
            % Perform the propagation again with higher number of discretization points, to better highligth the results:
            date0=[2028, 01, 01, 00, 00, 00];           % initial date for propagation
            t0=date2mjd2000(date0)*60*60*24;            % initial date for propagation [s]
            k_prop=500;                                 % number of orbits chosen for the propagation
            N_prop=100000;                              % chosen discretization for the propagation
            tspan=linspace(t0, t0+k_prop*T ,N_prop);    % timespan of the discretization
    
            % Propagation (Gauss's planetary equations):
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'RSW');
            dy_Gauss=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'RSW');
            [T_Gauss, Y_Gauss]=ode113(dy_Gauss, tspan, par_n, options);
            aGauss = Y_Gauss(:,1);
            eGauss = Y_Gauss(:,2);
            iGauss = Y_Gauss(:,3);
            OMGauss = Y_Gauss(:,4);
            wGauss = Y_Gauss(:,5);
            fGauss = Y_Gauss(:,6);
            
            % Filter the results coming from the propagation performed at point 3:
            T_short = T;                                                                % period of short-periodic oscillations [s]
            T_long = 15*T;                                                              % period of long-periodic oscillations [s]
            n_window_short = nearest(T_short/(sum(diff(T_Gauss)/(length(T_Gauss)-1)))); % number of points of the time window of the length of 1 T_short [s]
            n_window_long = nearest(T_long/(sum(diff(T_Gauss)/(length(T_Gauss)-1))));   % number of points of the time window of the length of 1 T_long [s]
            kep_filtered_long = movmean(Y_Gauss,n_window_short,1);                      % filtering of short-periodic oscillations (only secular and long periodic components remaining)
            kep_filtered_secular = movmean(Y_Gauss,n_window_long,1);                    % filtering of long-periodic oscillations (only secular components remaining)
            
            %% Plot filtered and unflitered keplerian elements:
            figure()
            subplot(3,1,1)
            plot ((T_Gauss-t0)/T, aGauss, 'Linewidth', 1e-6)
            hold on
            plot((T_Gauss-t0)/T, kep_filtered_long(:,1), 'Linewidth', 1e-6) 
            plot((T_Gauss-t0)/T, kep_filtered_secular(:,1), 'Linewidth', 1e-6,'Color','k') 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$a [km]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{Semi-major axis}', fontsize=12);
            grid on
            
            subplot(3,1,2)
            plot ((T_Gauss-t0)/T, eGauss, 'Linewidth', 1e-6)
            hold on
            plot((T_Gauss-t0)/T, kep_filtered_long(:,2), 'Linewidth', 1e-6) 
            plot((T_Gauss-t0)/T, kep_filtered_secular(:,2), 'Linewidth', 1e-6,'Color','k') 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$e [-]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{Eccentricity}', fontsize=12);
            grid on
            
            subplot(3,1,3)
            plot ((T_Gauss-t0)/T, rad2deg(iGauss), 'LineWidth',2)
            hold on
            plot((T_Gauss-t0)/T, rad2deg(kep_filtered_long(:,3)), 'LineWidth',1) 
            plot((T_Gauss-t0)/T, rad2deg(kep_filtered_secular(:,3)), 'Color','k', 'LineWidth',0.02) 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$i [deg]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{Inclination}', fontsize=12);
            grid on
            
            figure()
            subplot(3,1,1)
            plot ((T_Gauss-t0)/T, unwrap(rad2deg(OMGauss)), 'Linewidth', 1e-6)
            hold on
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_long(:,4))), 'Linewidth', 1e-6) 
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_secular(:,4))), 'Linewidth', 1e-6,'Color','k') 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\Omega [deg]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{Right Ascension of the Ascending Node}', fontsize=12);
            grid on
            
            subplot(3,1,2)
            plot ((T_Gauss-t0)/T, unwrap(rad2deg(wGauss)), 'Linewidth', 1e-6)
            hold on
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_long(:,5))), 'Linewidth', 1e-6) 
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_secular(:,5))), 'Linewidth', 1e-6,'Color','k') 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\omega [deg]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{Argument of Perigee}', fontsize=12);
            grid on
            
            subplot (3,1,3)
            plot ((T_Gauss-t0)/T, unwrap(rad2deg(fGauss)), 'Linewidth', 2)
            hold on
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_long(:,6))), 'Linewidth', 1) 
            plot((T_Gauss-t0)/T, unwrap(rad2deg(kep_filtered_secular(:,6))), 'Linewidth', 1e-6,'Color','k') 
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\theta [deg]$', fontsize=12);
            legend('Unfiltered','Long period evolution', 'Secular evolution');
            title('\bf{True Anomaly}', fontsize=12);
            grid on
                      
            % Algorithm to choose another task:
            choice5=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice5==2
                Perform_task='false';
            end 
        case 6
            %% Comparison with real data of a real celestial object:
            % Satellite used:
            % Name: SKYNET 1   NORAD Catalogue ID: 4250
            
            % Data taken from the ephemerides:
            options_eph = detectImportOptions('SKYNET_1.txt');
            options_eph.Delimiter=',';
            options_eph.VariableNamingRule = 'preserve';
            options_eph.VariableNames = {'JDTDB', 'Calendar Date (TDB)', 'EC', 'QR', 'IN', 'OM', 'W', 'Tp', 'N', 'MA', 'TA', 'A', 'AD', 'PR'};
            options_eph = setvartype(options_eph, 'Calendar Date (TDB)', 'string');
            ephemerides=readtable('SKYNET_1.txt', options_eph);
            a_sat_real = str2double(ephemerides{1:17545, 'A'});
            e_sat_real = str2double(ephemerides{1:17545, 'EC'});
            i_sat_real = deg2rad(str2double(ephemerides{1:17545, 'IN'}));
            OM_sat_real = deg2rad(str2double(ephemerides{1:17545, 'OM'}));
            w_sat_real = deg2rad(str2double(ephemerides{1:17545, 'W'}));
            f_sat_real = deg2rad(str2double(ephemerides{1:17545, 'TA'}));
            kep_sat_real=[a_sat_real, e_sat_real, i_sat_real, OM_sat_real, w_sat_real, f_sat_real];
            
            %The initial data for the propagation can be taken from ephemerides:
            a0_sat=a_sat_real(1);                                         % semi-major axis of the satellite at initial time
            e0_sat=e_sat_real(1);                                         % eccetricity of the satellite at initial time
            i0_sat=i_sat_real(1);                                         % inclination of the satellite at initial time
            OM0_sat=OM_sat_real(1);                                       % RAAN of the satellite at initial time
            w0_sat=w_sat_real(1);                                         % argument of pericenter of the satellite at initial time
            f0_sat=f_sat_real(1);                                         % true anomaly of the satellite at initial time
            kep0_sat = [a0_sat, e0_sat, i0_sat, OM0_sat, w0_sat, f0_sat]; % initial keplerian parameters

            % Define the time windows and spans:
            T_sat=2*pi*sqrt(a0_sat^3/mu_E);                               % period of the satellite
            initial_date = [2000, 02, 28, 00, 00, 00];                    % initial date for timespan
            t_in = date2mjd2000(initial_date);                            % initial date for timespan in MJD2000 (days)
            t_in_s = t_in*24*60*60;                                       % initial date for timespan in MJD2000 (seconds)
            final_date = [2002, 02, 28, 00, 00, 00];                      % final date for timespan
            t_fin = date2mjd2000(final_date);                             % final date for timespan in MJD2000 (days)
            t_fin_s = t_fin*24*60*60;                                     % final date for timespan in MJD2000 (seconds)
            step = 60*60;                                                 % time step in seconds for the propagation (chosen equal to the one of ephemerides)
            t_span_sat = t_in_s:step:t_fin_s;                             % time span for the propagation
            
            % Propagation of the perturbed model via Gauss's planetary equations:
            options = odeset('RelTol',1e-13,'AbsTol',1e-14);
            acc_pert=@(t, y) acc_pert_fun_J2_moon(t, y, mu_E, J2, mu_M, R_E, 'RSW');
            dy_Gauss=@(t,y) eq_motion(t, y, acc_pert(t,y), mu_E, 'RSW');
            [Time_sat, kep_sat]=ode113(dy_Gauss, t_span_sat, kep0_sat, options);
            
            % Perform the unwrap of the angular keplerian parameters and the conversion
            % to degree, this choice is made since there are no additional computation 
            % to do and it's useflu for plots:
            kep_sat_real(:,3:6)=rad2deg(unwrap(kep_sat_real(:,3:6)));
            kep_sat(:,3:6)=rad2deg(unwrap(kep_sat(:,3:6)));
            
            % Errors between propagation and real data from ephemerides:
            err_a_plot1 = abs(kep_sat(:,1)-kep_sat_real(:,1))./abs(kep0_sat(1));         % relative error of semi-major axis
            err_e_plot1 = abs(kep_sat(:,2)-kep_sat_real(:,2));                           % absolute error of eccentricity
            err_i_plot1 = abs(kep_sat(:,3)-kep_sat_real(:,3))/(360);                     % relative error of inclination
            err_OM_plot1 = abs(kep_sat(:,4)-kep_sat_real(:,4))/(360);                    % relative error of right ascension of ascending node
            err_om_plot1 = abs(kep_sat(:,5)-kep_sat_real(:,5))/(360);                    % relative error of argoument of periapsis
            err_f_plot1 = abs(kep_sat(:,6)-kep_sat_real(:,6))./abs(kep_sat_real(:,6));   % relative error of true anomaly
            
            % Maximum errors:
            max_err_a = max(err_a_plot1);
            max_err_e = max(err_e_plot1);
            max_err_i = max(err_i_plot1);
            max_err_OM = max(err_OM_plot1);
            max_err_om = max(err_om_plot1);
            max_err_th = max(err_f_plot1);
            
            % Plot the comparison between the propagation method and the real data from
            % ephemerides and related error:
            figure()
            subplot(2,2,1)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,1),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,1),'LineWidth',0.5)
            title('\bf{Evolution of $a$ for Gauss propagation and ephemerides [km]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period',fontsize=12);
            ylabel('$a [km]$', fontsize=12);
            legend('Evolution of $a_{Gauss}$', 'Evolution of $a_{ephemerides}$')
            grid on
         
            subplot(2,2,2)
            semilogy((Time_sat-t_in_s)/T_sat,err_a_plot1,'LineWidth',0.5)
            title('\bf{Relative error of $a$}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|a_{Gauss} - a_{eph}|/|a_0|$ ', fontsize=14);
            grid on
  
            subplot(2,2,3)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,2),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,2),'LineWidth',0.5)
            title('\bf{Evolution of $e$ for Gauss propagation and ephemerides [-]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$e$ [-]', fontsize=12);
            legend('Evolution of $e_{Gauss}$', 'Evolution of $e_{ephemerides}$')
            grid on
            
            subplot(2,2,4)
            semilogy((Time_sat-t_in_s)/T_sat,err_e_plot1,'LineWidth',0.5)
            title('\bf{Absolute error of eccentricity}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|e_{Gauss} - e_{eph}|$', fontsize=14);
            grid on
       
            figure()
            subplot(2,2,1)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,3),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,3),'LineWidth',0.5)
            legend('Evolution of $i_{Gauss}$', 'Evolution of $i_{ephemerides}$')
            title('\bf{Evolution of $i$ for Gauss propagation and ephemerides [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$i$ [deg]', fontsize=12);
            grid on
           
            subplot(2,2,2)
            semilogy((Time_sat-t_in_s)/T_sat,err_i_plot1,'LineWidth',0.5)
            title('\bf{Relative error of inclination}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|i_{Gauss} - i_{eph}|/360^{\circ}$', fontsize=14);
            grid on
      
            subplot(2,2,3)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,4),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,4),'LineWidth',0.5)
            legend('Evolution of $\Omega_{Gauss}$', 'Evolution of $\Omega_{ephemerides}$')
            title('\bf{Evolution of $\Omega$ for Gauss propagation and ephemerides [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\Omega [deg]$', fontsize=12);
            grid on
            
            subplot(2,2,4)
            semilogy((Time_sat-t_in_s)/T_sat,err_OM_plot1,'LineWidth',0.5)
            title('\bf{Relative error of the right ascension of the ascending node}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\Omega_{Gauss} - \Omega_{eph}|/360^{\circ}$', fontsize=14);
            grid on

            figure()
            subplot(2,2,1)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,5),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,5),'LineWidth',0.5)
            legend('Evolution of $\omega_{Gauss}$', 'Evolution of $\omega_{ephemerides}$')
            title('\bf{Evolution of $\omega$ for Gauss propagation and ephemerides [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\omega [deg]$', fontsize=12);
            grid on

            subplot(2,2,2)
            semilogy((Time_sat-t_in_s)/T_sat,err_om_plot1,'LineWidth',0.5)
            title('\bf{Relative error of the argument of periapsis}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\omega_{car} - \omega_{Ephemerides}|/360^{\circ}$ ', fontsize=14);
            grid on
            
            subplot(2,2,3)
            plot((Time_sat-t_in_s)/T_sat,kep_sat(:,6),'LineWidth',0.5)
            hold on
            plot((Time_sat-t_in_s)/T_sat,kep_sat_real(:,6),'LineWidth',0.5)
            legend('Evolution of $\theta_{Gauss}$', 'Evolution of $\theta_{ephemerides}$')
            title('\bf{Evolution of $\theta$ Gauss propagation and ephemerides [deg]}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$\theta$ [deg]', fontsize=12);
            grid on

            subplot(2,2,4)
            semilogy((Time_sat-t_in_s)/T_sat,err_f_plot1,'LineWidth',0.5)
            title('\bf{Relative error of the true anomaly}', fontsize=12);
            xlabel('$N^{\circ}$ of orbit period', fontsize=12);
            ylabel('$|\theta_{Gauss} - \theta_{eph}|/|\theta_{eph}|$ ', fontsize=14);
            grid on

            % Algorithm to choose another task:
            choice6=input(['Do you want to perform another task? \n'...
                '1: Yes \n' ...
                '2: No \n']);
            if choice6==2
                Perform_task='false';
            end 
        case 7
            Perform_task='false';
        otherwise 
            disp('Invalid task identifier value')
    end 
end