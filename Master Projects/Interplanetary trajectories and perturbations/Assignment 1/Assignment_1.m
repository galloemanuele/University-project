%% Main script for Interplanetary Explorer Mission (Assignment 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       ORBITAL MECHANICS                                 %
%                    Academic year 2023/2024                              %
%                    M.Sc. Space Engineering                              %
%                     Politecnico di Milano                               %
%                                                                         %
%              Interplanetary Explorer Mission assignment                 %     
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

%% DATA
% Departure:          Mercury 
% Fly-By:             Earth 
% Arrival:            Asteroid N.21 
% Earliest Departure: 00:00:00 01/01/2028 
% Latest Arrival:     00:00:00 01/01/2058
% Path to the functions:

% Path to the functions and textures:
addpath(strcat(pwd,'/Functions'));
addpath(strcat(pwd,'/Textures'));
addpath(strcat(pwd,'/Given functions'));
addpath(strcat(pwd,'/Given functions/timeConversion'));
clear
clc
close 

mu_Sun = astroConstants(4);
mu_E=astroConstants(13);
% Definition of the departure and arrival planets
dep_planet = 1; % Mercury
flyby_planet = 3; % Earth
arr_body = 21; % Asteroid N.21

% Definition of Launch and Arrival Windows & Time Conversion
dep_min = [2028, 01, 01, 00, 00, 00]; % January 1st, 2028 h 00:00:00
dep_min = date2mjd2000(dep_min); 

arr_max = [2058, 01, 01, 00, 00, 00]; % January 1st, 2058 h 00:00:00
arr_max = date2mjd2000(arr_max);

% Calculation of Planets Periods
[kep_dep, ~] = uplanet(0, dep_planet);
T_dep = 2*pi*sqrt(kep_dep(1)^3/mu_Sun)*1/86400; % Orbital period of departure planet

[kep_flyby, ~] = uplanet(0, flyby_planet);
T_Flyby = 2*pi*sqrt(kep_flyby(1)^3/mu_Sun)*1/86400; % Orbital period of Fly-By planet

[kep_arr, ~] = ephNEO(0, arr_body);
T_arr = 2*pi*sqrt(kep_arr(1)^3/mu_Sun)*1/86400; % Orbital period of arrival body

%Hohmann like tranfer
[Dvi1,Dvf1,~,Dt1,~,~]= bitangent(kep_dep(1),kep_dep(2),0,kep_flyby(1),kep_flyby(2),0,'pa',mu_Sun);
[Dvi2,Dvf2,~,Dt2,~,~] = bitangent(kep_flyby(1),kep_flyby(2),0,kep_arr(1),kep_arr(2),0,'pa',mu_Sun); 
dt1=Dt1*1/86400;
dt2=Dt2*1/86400;

%% ToFs estimation process
fprintf("To use this program it is necessary to have downloaded the Optimization Toolbox \nand the Statistics and Machine Learning ToolBox")
var=input("\nChoose \n 1 Visualize Porkchop plot with fixed flyby date and the comparison with the localization of the minima through the grid" +...
    " \n 2 Run the process of localization of the local minima to visualize the process of ToFs determination "  +...
    "\n 3 Run directly in the optimisation process and in the showing of the results \n " + ...
    "4 Exit the program \n");
switch var
    case 1
    t2=date2mjd2000([2035,2 ,1, 00, 00, 00]);
    FixedFB_Porkchop(t2,dep_min,arr_max,T_arr,T_dep)
    case 2
    t2=linspace(dep_min,arr_max,200);
    P=BuildGrid(t2,dep_min,arr_max,T_dep,T_arr);
    
    for i=1:length(t2) %Fixed flyby date
        t1=P(:,1,i);
        t3=P(:,2,i); 
            for j=1:length(t1)
                if not(isnan(t1(j))) && not(isnan(t3(j)))
                    Dvtot_mat3D(j,j,i) = dV_interplanetary_unc(t1(j),t2(i),t3(j));
           
            %For every flyby date i I have a matrix of cost depending on
            %departure and arrival.
                else
                    Dvtot_mat3D(j,j,i)=NaN;
                end
            end
    end
    
    [m,n,l]=size(Dvtot_mat3D);
    for i=1:l
        for j=1:m
            for k=1:n
                if  Dvtot_mat3D(j,k,i)==0
                    Dvtot_mat3D(j,k,i)=NaN;
                end
            end
        end
    end
    t_flyby=zeros(length(t2),3);
    Dv_fbtest=zeros(length(t2),1);
    
    for i=1:l
      t1=P(:,1,i);
      t3=P(:,2,i); 
      [Dv_fbtest(i),k_opt]=min(min(Dvtot_mat3D(:,:,i)));
      [~,j_opt]=min(Dvtot_mat3D(:,k_opt,i));
      if isnan(Dv_fbtest(i,:)) 
         t_flyby(i,:)=[0 0 0];
        
      else
      t_flyby(i,:)=[t1(j_opt) t2(i) t3(k_opt)];
      end
    
    end
    ToF1=zeros(length(t2),1);
    ToF2=zeros(length(t2),1);
    ToF_tot=zeros(length(t2),1);
    for i=1:l
        
        ToF1(i)=t_flyby(i,2)-t_flyby(i,1);
        ToF2(i)=t_flyby(i,3)-t_flyby(i,2);
        ToF_tot(i)=ToF2(i)+ToF1(i);
    end
    
    %% Calculate statistics for ToFs
    ToF1(ToF1==0)=NaN;
    ToF2(ToF2==0)=NaN;
    
    Dv_min_unc=min(min(Dv_fbtest));
    ToF1_mean = mean(ToF1, 'omitnan');
    ToF1_min = ToF1_mean - std(ToF1, 'omitnan');
    ToF1_max = ToF1_mean + std(ToF1, 'omitnan');
    
    ToF2_mean = mean(ToF2, 'omitnan');
    ToF2_min = ToF2_mean - std(ToF2, 'omitnan');
    ToF2_max = ToF2_mean + std(ToF2, 'omitnan');
    
    
    % Generate data for kernel density estimation
    x_values = linspace(min([min(ToF1), min(ToF2)]), max([max(ToF1), max(ToF2)]), 1000);
    
    % Kernel density estimation for ToF1
    pdf1 = ksdensity(ToF1, x_values);
    
    % Kernel density estimation for ToF2
    pdf2 = ksdensity(ToF2, x_values);
    
    % Plot histograms and kernel density estimates in subplots
    figure;
    
    % Subplot 1 for ToF1
    subplot(2, 1, 1);
    histogram(ToF1, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'w', 'DisplayName', 'ToF1 Histogram');
    hold on;
    plot(x_values, pdf1, 'r', 'LineWidth', 2, 'DisplayName', 'ToF1 KDE');
    line([ToF1_mean, ToF1_mean], [0, max(pdf1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Mean ToF1');
    line([ToF1_min, ToF1_min], [0, max(pdf1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Min ToF1');
    line([ToF1_max, ToF1_max], [0, max(pdf1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Max ToF1');
    title('ToF1 Distribution');
    xlabel('ToF1 Values');
    ylabel('Probability Density');
    legend('show');
    hold off;
    
    % Subplot 2 for ToF2
    subplot(2, 1, 2);
    histogram(ToF2, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'w', 'DisplayName', 'ToF2 Histogram');
    hold on;
    plot(x_values, pdf2, 'b', 'LineWidth', 2, 'DisplayName', 'ToF2 KDE');
    line([ToF2_mean, ToF2_mean], [0, max(pdf2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Mean ToF2');
    line([ToF2_min, ToF2_min], [0, max(pdf2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Min ToF2');
    line([ToF2_max, ToF2_max], [0, max(pdf2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Max ToF2');
    title('ToF2 Distribution');
    xlabel('ToF2 Values');
    ylabel('Probability Density');
    legend('show');
    hold off;
    
    sgtitle('Histograms and Kernel Density Estimates of ToF1 and ToF2');
    case 3
    %% Genetic algorithm
    var=input("Choose \n 1  Run the genetic algorithm with the ToFs estimated by using the local minima individuation process" +...
    " \n 2  Run the genetic algorithm with the ToFs estimated by using the Hohmann like transfer \n");
    %From the observation of the graphs
    if var==1
        %Values estimated previously
        ToF1_ga_min=2.337876744449154e+02;
        ToF1_ga_max=4.999712649956725e+02;
        ToF2_ga_min=2.110408565949621e+02;
        ToF2_ga_max=1.337507621287187e+03;
      
        ToF1=[ToF1_ga_min;ToF1_ga_max];
        ToF2=[ToF2_ga_min;ToF2_ga_max];
    else
        ToF1=[dt1*0.1;dt1*2];
        ToF2=[dt1*0.1;dt2*2];
    end

    launch_w=[dep_min;arr_max];
    
    
    options_ga = optimoptions('ga','MaxGeneration',100,'PopulationSize',1000,...
        'FunctionTolerance',0.01,...
        'PlotFcn',"gaplotbestf",'Display','off');
    rng default                               % For riproducibility
    N_runs = 5;                               % Number of runs
    lb_ga = [launch_w(1) ToF1(1) ToF2(1)];              % Lower boundary   
    ub_ga = [launch_w(end) ToF1(end) ToF2(end)];        % Upper boundary
    
    % Initialization of output array
    Dv_ga_vect = ones(N_runs,1);              % minimum dv array obtained via Genetic Algorithm
    x_ga = zeros(N_runs,3);                   % minimum dv time array obtained via Genetic Algorithm
    
    fprintf('Genetic Algorithm running...\n');
    tic
    for i = 1:N_runs
    
        [x_ga(i,:),Dv_ga_vect(i)] = ga(@(x) dV_interplanetary_ToF(x(1),x(2),x(3)),3,[],[],[],[],lb_ga,ub_ga,[],options_ga);
    end
    fprintf('\nFinished. \n');
    toc
    
    [dv_ga,ii] = min(Dv_ga_vect);
    t_ga =[x_ga(ii,1) x_ga(ii,1)+x_ga(ii,2) x_ga(ii,1)+x_ga(ii,2)+x_ga(ii,3)] ;  
    % time array of the minimum dv obtained via Genetic Algorithm
    ToF1_ga=x_ga(ii,2);
    ToF2_ga=x_ga(ii,3);
    
    fprintf('Results:\n');
    fprintf('dv_ga = %.4f km/s\n', dv_ga)
    fprintf('dt_ga = %.4f days\n', ToF1_ga+ToF2_ga)
    fprintf('Departure date: %s\n', date2string(mjd20002date(t_ga(1))))
    fprintf('Flyby Date %s\n', date2string(mjd20002date(t_ga(2))))
    fprintf('Arrival date: %s\n', date2string(mjd20002date(t_ga(3))))
    fprintf('----------------------------------------------\n');
    
    %% Fmincon
    t1_fmin=linspace(t_ga(1)-1,t_ga(1)+1,2);                                              
    lb_fmincon = [t1_fmin(1) ToF1_ga-5 ToF2_ga-5];           % lower boundary
    ub_fmincon = [t1_fmin(end) ToF1_ga+5 ToF2_ga+5];         % upper boundary
    x0 = x_ga(ii,:);                                         % starting point
    options = optimoptions('fmincon','Algorithm','sqp','PlotFcn',"optimplotfval",...
        'FunctionTolerance',1e-10,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,...
        'display','iter');
    
    fprintf('Gradient-based optimization running...\n');
    
    tic
    [t_fmincon,Dv_fmincon] = fmincon(@(x) dV_interplanetary_ToF(x(1),x(2),x(3)),x0,[],[],[],[],lb_fmincon,ub_fmincon,[],options);
    fprintf('Finished. \n')
    toc
    
    %Dates
    t_fmincon=[t_fmincon(1)  t_fmincon(1)+t_fmincon(2)  t_fmincon(1)+t_fmincon(2)+t_fmincon(3)];
    fprintf('Results:\n');
    fprintf('Dv_fmincon = %.4f km/s\n', Dv_fmincon)
    fprintf('dt_fmincon = %.4f days\n', t_fmincon(3)-t_fmincon(1))
    fprintf('Departure date: %s\n', date2string(mjd20002date(t_fmincon(1))))
    fprintf('Flyby Date %s\n', date2string(mjd20002date(t_fmincon(2))))
    fprintf('Arrival date: %s\n', date2string(mjd20002date(t_fmincon(3))))
    fprintf('----------------------------------------------\n');
    
    %% Grid search
    t1_grid=linspace(t_ga(1)-1,t_ga(1)+2,100);
    ToF1_grid=linspace(ToF1_ga-1,ToF1_ga+1,100);
    ToF2_grid=linspace(ToF2_ga-1,ToF2_ga+2,100);
    Dvtot_mat = zeros(length(t1_grid),length(ToF1_grid),length(ToF2_grid));
    fprintf('Grid search running...\n');
    for i = 1:length(t1_grid)
        fprintf("\n iteration %.1d ",i)
        for j = 1:length(ToF1_grid)
            for k = 1:length(ToF2_grid)
    
            Dvtot_mat(i,j,k) = dV_interplanetary_ToF(t1_grid(i),ToF1_grid(j),ToF2_grid(k));
    
            end
        end
    end 
    
    [Dv_grid,loc] = min(Dvtot_mat(:));
    [i_opt,j_opt,k_opt] = ind2sub(size(Dvtot_mat),loc);
    %Dates
    t_grid = [t1_grid(i_opt) t1_grid(i_opt)+ToF1_grid(j_opt) t1_grid(i_opt)+ToF1_grid(j_opt)+ToF2_grid(k_opt)]; 
    
    fprintf('Results:\n');
    fprintf('dv_Grid = %.4f km/s\n', Dv_grid)
    fprintf('dt_Grid = %.4f days\n', ToF1_grid(j_opt)+ToF2_grid(k_opt))
    fprintf('Departure date: %s\n', date2string(mjd20002date(t_grid(1))))
    fprintf('Flyby Date %s\n', date2string(mjd20002date(t_grid(2))))
    fprintf('Arrival date: %s\n', date2string(mjd20002date(t_grid(3))))
    fprintf('----------------------------------------------\n');
    %% Results
    Results_Dv=[Dv_fmincon,Dv_grid];
    Results_Dt=[t_fmincon;t_grid];
    [~,index]=min(Results_Dv);
    t_opt=Results_Dt(index,:);
    ToF1_opt=t_opt(2)-t_opt(1);
    ToF2_opt=t_opt(3)-t_opt(2);
    
    
    [Dv_opt,r1_opt,r2_opt,r3_opt,Dv1_vect,Dv2_vect,~,vt1_in,vt2_in,vt1_out,vt2_out,v_inf_1,v_inf_2] = dV_interplanetary_ToF(t_opt(1),ToF1_opt,ToF2_opt);
    dV_flyby=norm(vt2_in-vt1_out);
    [Dv_pow,r_p,a_hyp1,e_hyp1,a_hyp2,e_hyp2,v_p1,v_p2]=PoweredFlyBy(v_inf_1,v_inf_2);
    %ToF in the SOI
    kep_fb_opt = uplanet(t_opt(2), flyby_planet);
    [rrE_flyby,vvE_flyby] =  kep2car(kep_fb_opt(1), kep_fb_opt(2), kep_fb_opt(3), kep_fb_opt(4), kep_fb_opt(5),kep_fb_opt(6), mu_Sun);
    ToF_SOI=ToF_flyby(a_hyp1,e_hyp1,a_hyp2,e_hyp2,rrE_flyby);
    % Lambert arcs characterization
    Kep_Lambert1_in = car2kep(r1_opt, vt1_in, mu_Sun);
    Kep_Lambert1_fin = car2kep(r2_opt, vt1_out, mu_Sun);
    Kep_Lambert2_in = car2kep(r2_opt, vt2_in, mu_Sun);
    Kep_Lambert2_fin = car2kep(r3_opt, vt2_out, mu_Sun);
    format long
    fprintf('Results:\n');
    fprintf('dv_opt = %.4f km/s\n', Dv_opt)
    fprintf('dv_pow = %.7f km/s\n', Dv_pow)
    fprintf('Departure date: %s\n', date2string(mjd20002date(t_opt(1))))
    fprintf('Flyby Date %s\n', date2string(mjd20002date(t_opt(2))))
    fprintf('Arrival date: %s\n', date2string(mjd20002date(t_opt(3))))
    fprintf('ToF in the SOI: %4f [h] \n', ToF_SOI/(60*60))
    fprintf('----------------------------------------------\n');

    %% Plots
    var=input("Choose \n 1 Non animated plots \n 2 Animated plot  \n 3 Exit the program\n");
    
    
    switch var
        case 1
        
        ToF1_s=ToF1_opt*24*3600;
        ToF2_s=ToF2_opt*24*3600;
        % Set time span
        dt = 100; % Step size [s]
        AU=astroConstants(2);
        figure(1)
        % Trajectory of Mercury
        tspan= 0:dt:ToF1_s;
        kep0 = uplanet(t_opt(1), dep_planet);
        [rr0_dep, vv0_dep] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_Sun);
        y0 = [rr0_dep; vv0_dep];
        [~,Y_dep] = orbitPropagator(y0,tspan,"up",mu_Sun);
        Traj1=plot3( Y_dep(:,1)./AU, Y_dep(:,2)./AU, Y_dep(:,3)./AU, '-','LineWidth',1,'Color','b');
        hold on
        
        % Complete initial orbit
        T = 2*pi*sqrt(kep0(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [T_depc,Y_depc]=orbitPropagator(y0,tspan,"up",mu_Sun);
        Initial=plot3( Y_depc(:,1)./AU, Y_depc(:,2)./AU, Y_depc(:,3)./AU, '--','LineWidth',1,'Color','b');
        
        % Trajectory of Earth
        tspan = 0:dt:ToF1_s;
        kep_flyby = uplanet(t_opt(1), flyby_planet);
        [rr0_fb, vv0_fb] = kep2car(kep_flyby(1), kep_flyby(2), kep_flyby(3), kep_flyby(4), kep_flyby(5), kep_flyby(6), mu_Sun);
        y0 = [rr0_fb; vv0_fb];
        [~,Y_FB]=orbitPropagator(y0,tspan,"up",mu_Sun);
        Traj2=plot3( Y_FB(:,1)./AU, Y_FB(:,2)./AU, Y_FB(:,3)./AU, '-','LineWidth',1,'Color','red');
        
        % Complete final orbit
        T = 2*pi*sqrt(kep_flyby(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [~, Y_FBC] = orbitPropagator(y0,tspan,"up",mu_Sun);
        Fb=plot3( Y_FBC(:,1)./AU, Y_FBC(:,2)./AU, Y_FBC(:,3)./AU, '--','LineWidth',1,'Color','red');
        
        % Trajectory of Asteroid N.21
        tspan = 0:dt:ToF1_s+ToF2_s;
        kep_arr = ephNEO(t_opt(1), arr_body);
        [rr0_arr, vv0_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_Sun);
        y0 = [rr0_arr; vv0_arr];
        [~,Y_arr]=orbitPropagator(y0,tspan,"up",mu_Sun);
        Traj3=plot3( Y_arr(:,1)./AU, Y_arr(:,2)./AU, Y_arr(:,3)./AU, '-','LineWidth',1,'Color','green');
        
        % Complete flyby orbit
        T = 2*pi*sqrt(kep_arr(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [~, Y_arrC] = orbitPropagator(y0,tspan,"up",mu_Sun);
        Final=plot3( Y_arrC(:,1)./AU, Y_arrC(:,2)./AU, Y_arrC(:,3)./AU, '--','LineWidth',1,'Color','green');
        
        % Lambert arc 1
        tspan = 0:dt:ToF1_s;
        y0 = [r1_opt; vt1_in];
        [~, Yt1] = orbitPropagator(y0,tspan,"up",mu_Sun);
        Lambert1=plot3( Yt1(:,1)./AU, Yt1(:,2)./AU, Yt1(:,3)./AU, '-','LineWidth',2,'Color','yellow');
        
        % Lambert arc 2
        tspan = 0:dt:ToF2_s;
        y0 = [r2_opt; vt2_in];
        [~, Yt2] =  orbitPropagator(y0,tspan,"up",mu_Sun);
        Lambert2=plot3( Yt2(:,1)./AU, Yt2(:,2)./AU, Yt2(:,3)./AU, '-','LineWidth',2,'Color','cyan');
        
        %Sun
        Plot_Planet(10,1.4378e-07);
        %Mercury
        Plot_Planet(1,1.2791e-05,r1_opt(1)/AU,r1_opt(2)/AU,r1_opt(3)/AU);
        hold on
        %Earth
        Plot_Planet(3,1.0557e-05,r2_opt(1)/AU,r2_opt(2)/AU,r2_opt(3)/AU);
        %Asteroid
        Plot_Planet(11,1.1791e-05,r3_opt(1)/AU,r3_opt(2)/AU,r3_opt(3)/AU);
        grid on

        title('Interplanetary transfer strategy');
        xlabel('X [AU]');
        ylabel('Y [AU]');
        zlabel('Z [AU]');
        legend( ...
            'Mercury Initial Orbit', 'Mercury Trajectory','Earth Trajectory', 'Earth Flyby Orbit', ...
            'Asteroid Trajectory', 'Asteroid Final Orbit', 'Lambert Arc 1', 'Lambert Arc 2','Sun','Mercury','Earth','Asteroid 21');
          
        legend('Location', 'Best');
        legend('boxoff');
       
        %% Flyby plot------------------------------------------------------
        figure(2)
        R_E=astroConstants(23);
        %Normal to the manoeuvre plane
        u = cross(v_inf_1,v_inf_2)/norm(cross(v_inf_1,v_inf_2));
        
        %Planetocentric arcs in Earth centred frame parallel to HECI
        del_minus=asin(1/e_hyp1);
        figure(2)
        Earth3D(1)
        hold on; grid on;
        rotation_angle = del_minus + pi/2;
        a_vector =-rodrigues(v_inf_1, u, rotation_angle); %with this operation you find the opposite vector, so you put a -
        r_p_Point = r_p * a_vector/norm(a_vector);
        b_vector= rodrigues(v_inf_1, u, del_minus);
        b_vector= b_vector/norm(b_vector);
        v_p1_vect=b_vector*v_p1;
        v_p2_vect=b_vector*v_p2;
        
        hold on; grid on
        y_b=[r_p_Point;v_p1_vect];
        y_a=[r_p_Point;v_p2_vect];
        tspan=linspace(0,-1000,1000);
        [~,Y_b]=orbitPropagator(y_b,tspan,"up",mu_E);
        IncomingHyperbola=plot3(Y_b(:,1)./R_E,Y_b(:,2)./R_E,Y_b(:,3)./R_E,'LineWidth',1,'Color','blue');
        
        tspan=-tspan;
        [~,Y_a]=orbitPropagator(y_a,tspan,"up",mu_E);
        OutgoingHyperbola=plot3(Y_a(:,1)./R_E,Y_a(:,2)./R_E,Y_a(:,3)./R_E,'LineWidth',1,'Color','red');
        % scatter3(Y_a(1,1)./R_E,Y_a(1,2)./R_E,Y_a(1,3)./R_E, 'green','filled','o')
        
        legend([IncomingHyperbola, OutgoingHyperbola], 'Incoming hyperbola', 'Outgoing hyperbola');
        title('Fly-by');
        xlabel('$x \, [R_e]$', 'Interpreter', 'latex');
        ylabel('$y \, [R_e]$', 'Interpreter', 'latex');
        zlabel('$z \, [R_e]$', 'Interpreter', 'latex');

        
        
        case 2
      
        ToF1_s=ToF1_opt*24*3600;
        ToF2_s=ToF2_opt*24*3600;
        % Set time span
        dt = 100; % Step size [s]
        AU=astroConstants(2);
        hFig=figure(1);
        set(hFig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        % Trajectory of Mercury
        tspan= 0:dt:ToF1_s;
        kep0 = uplanet(t_opt(1), dep_planet);
        [rr0_dep, vv0_dep] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_Sun);
        y0 = [rr0_dep; vv0_dep];
        [~,Y_dep] = orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_dep(:,1)./AU, Y_dep(:,2)./AU, Y_dep(:,3)./AU, '-','LineWidth',1,'Color','b');
        hold on
        
        % Complete initial orbit
        T = 2*pi*sqrt(kep0(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [~,Y_depc]=orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_depc(:,1)./AU, Y_depc(:,2)./AU, Y_depc(:,3)./AU, '--','LineWidth',1,'Color','b');
        
        % Trajectory of Earth
        tspan = 0:dt:ToF1_s;
        kep_flyby = uplanet(t_opt(1), flyby_planet);
        [rr0_fb, vv0_fb] = kep2car(kep_flyby(1), kep_flyby(2), kep_flyby(3), kep_flyby(4), kep_flyby(5), kep_flyby(6), mu_Sun);
        y0 = [rr0_fb; vv0_fb];
        [~,Y_FB]=orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_FB(:,1)./AU, Y_FB(:,2)./AU, Y_FB(:,3)./AU, '-','LineWidth',1,'Color','red');
        
        % Complete final orbit
        T = 2*pi*sqrt(kep_flyby(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [~, Y_FBC] = orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_FBC(:,1)./AU, Y_FBC(:,2)./AU, Y_FBC(:,3)./AU, '--','LineWidth',1,'Color','red');
        
        % Trajectory of Asteroid N.21
        tspan = 0:dt:(ToF1_s+ToF2_s);
        kep_arr = ephNEO(t_opt(1), arr_body);
        [rr0_arr, vv0_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_Sun);
        y0 = [rr0_arr; vv0_arr];
        [~,Y_arr]=orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_arr(:,1)./AU, Y_arr(:,2)./AU, Y_arr(:,3)./AU, '-','LineWidth',1,'Color','green');
        
        % Complete flyby orbit
        T = 2*pi*sqrt(kep_arr(1)^3/mu_Sun);
        tspan = 0:dt:T;
        [~, Y_arrC] = orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Y_arrC(:,1)./AU, Y_arrC(:,2)./AU, Y_arrC(:,3)./AU, '--','LineWidth',1,'Color','green');
        
        % Lambert arc 1
        tspan = 0:dt:ToF1_s;
        y0 = [r1_opt; vt1_in];
        [~, Yt1] = orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Yt1(:,1)./AU, Yt1(:,2)./AU, Yt1(:,3)./AU, '-','LineWidth',2,'Color','yellow');
        
        % Lambert arc 2
        tspan = 0:dt:ToF2_s;
        y0 = [r2_opt; vt2_in];
        [~, Yt2] =  orbitPropagator(y0,tspan,"up",mu_Sun);
        plot3( Yt2(:,1)./AU, Yt2(:,2)./AU, Yt2(:,3)./AU, '-','LineWidth',2,'Color','cyan');
        
        title('Interplanetary transfer strategy');
        xlabel('X [AU]');
        ylabel('Y [AU]');
        zlabel('Z [AU]');

        %Sun
        Plot_Planet(10,1.4378e-07);
        %Mercury
        Plot_Planet(1,1.2791e-05,r1_opt(1)/AU,r1_opt(2)/AU,r1_opt(3)/AU);
        hold on
        %Earth
        Plot_Planet(3,1.0557e-05,r2_opt(1)/AU,r2_opt(2)/AU,r2_opt(3)/AU);
        %Asteroid
        Plot_Planet(11,1.1791e-05,r3_opt(1)/AU,r3_opt(2)/AU,r3_opt(3)/AU);

       legend(  'Mercury Initial Orbit', 'Mercury Trajectory','Earth Trajectory', 'Earth Flyby Orbit', ...
        'Asteroid Trajectory', 'Asteroid Final Orbit', 'Lambert Arc 1', 'Lambert Arc 2','Sun','Mercury','Earth','Asteroid 21');
        legend('Location', 'Best');
        legend('boxoff');
        
        grid on
        % Plot initial position of the satellite for Lambert arc 1
        satellitePosition1 = plot3(Yt1(1,1)./AU, Yt1(1,2)./AU, Yt1(1,3)./AU, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Satellite - Lambert 1');
        
        % Plot initial position of the satellite for Lambert arc 2
        satellitePosition2 = plot3(Yt2(1,1)./AU, Yt2(1,2)./AU, Yt2(1,3)./AU, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Satellite - Lambert 2');
        
        
        
        % Add a legend
        legend('Location', 'Best');
        legend('boxoff');
        
        % Add time annotation
        timeAnnotation = annotation('textbox', [0.75, 0.85, 0.1, 0.1], 'String', '', 'EdgeColor', 'none', 'FontSize', 12);
        
        % Set axis limits
        axis equal;
        axis tight;

        
        tspan1=linspace(0,ToF1_opt,100);
        tspan2=linspace(ToF1_opt,ToF2_opt,100);
        [~, Yt1] =  orbitPropagator([r1_opt; vt1_in],linspace(0,ToF1_s,100),"up",mu_Sun);
        [~, Yt2] =  orbitPropagator([r2_opt; vt2_in],linspace(0,ToF2_s,100),"up",mu_Sun);


        % videoFile = 'interplanetary_transfer_animation.mp4';
        % writerObj = VideoWriter(videoFile, 'MPEG-4');
        % open(writerObj);
        % Loop through time steps and update satellite positions and time annotation
        for i = 1:length(tspan1)
            % Update satellite position for Lambert arc 1
            set(satellitePosition1, 'XData', Yt1(i,1)./AU, 'YData', Yt1(i,2)./AU, 'ZData', Yt1(i,3)./AU);
        
            % Update time annotation
            current_time = tspan1(i); % Update with the correct time variable for Lambert arc 1
            set(timeAnnotation, 'String', sprintf('Time - Lambert 1: %.2f days', current_time));
        
            % Redraw the plot
            drawnow;
            
             % Write each frame to the video
            % writeVideo(writerObj, getframe(gcf));
            % Pause to control the animation speed
            pause(1e-20);
        end
        
        % Loop through time steps and update satellite position for Lambert arc 2
        for i = 1:length(Yt2)
            % Update satellite position for Lambert arc 2
            set(satellitePosition2, 'XData', Yt2(i,1)./AU, 'YData', Yt2(i,2)./AU, 'ZData', Yt2(i,3)./AU);
        
            % Update time annotation
            current_time = tspan2(i)+tspan1(end); % Update with the correct time variable for Lambert arc 2
            set(timeAnnotation, 'String', sprintf('Time - Lambert 2: %.2f days', current_time));
        
            % Redraw the plot
            drawnow;
            % Write each frame to the video
            % writeVideo(writerObj, getframe(gcf));
            % Pause to control the animation speed
            pause(1e-20);
        end
       
     % close(writerObj);
        case 3
    end
     case 4

end
