% Blowdown versione 5:

clear
close all
clc
%%
% Load design results:
load("RISULTATI_ITERAZIONE0.mat")

% The iteration i is referred to the instant in the combustion chamber,
% which is useful to calculate the performances.

% Initialize the cycle:
i=1;
t(1)=0; % [s]
% dt = 1; 
tmax=3600; % [s]     % Aim: dt must be small
options = optimoptions('fsolve', 'Display', 'off');


% Reallocate data coming from RISULTATI_ITERAZIONE0.mat: % hp liquido incomprimibile
blowdown.Pc(i)=Pc_in;
blowdown.Pc_f=Pc_in;
blowdown.M_esp_ox(i)=0;
% blowdown.M_esp_ox(i+1)=0;
blowdown.M_esp_f(i)=0;
% blowdown.M_esp_f(i+1)=0;
blowdown.r_of(i)=o_f;
blowdown.ve(i)=ve_in;
blowdown.m_ox(i)=m_OX;
blowdown.m_f(i)=m_FU;
blowdown.v_ox(i)=15; % m/s
blowdown.v_f(i)=15; % m/s
blowdown.V_Heox(i)=Vi_He_ox; % initial volume of helium in the oxigen tank
blowdown.V_Hef(i)=Vi_He_fu; % initial volume of helium in the fuel tank
blowdown.P_Heox(i)=P_tank_ox_in;
blowdown.P_Hef(i)=P_tank_f_in;
blowdown.V_ox(i)=Vi_ox;
blowdown.V_f(i)=Vi_fu;
blowdown.Mox(i)=Mi_ox; % vettore massa residua oxigen
blowdown.Mfu(i)=Mi_fu;
blowdown.T_He_ox(i)=Ti_He_ox;% K
blowdown.T_He_f(i)=Ti_He_fu; % K
blowdown.M_esp_tot_f=0;
blowdown.M_esp_tot_ox=0;
blowdown.mp(i)=mp_in;
blowdown.x(i)=0.15;
blowdown.y(i)=0.15;

performance.Itot=0;
performance.gamma_gc(i)=gamma_in;
performance.Mm(i)=Mm_in;
performance.Tf(i)=Tf_in;          % T flame 
performance.c_star(i)= c_star_in;
performance.Isp(i)=I_sp;
performance.cf(i)=c_f;
performance.M_cc(i)=M_cc_in;
performance.T_2D(i)=T_in_real;
performance.T(i)=T_in_ideal;
performance.rho(i)=rho_in;
performance.v_cc(i)=v_cc_in;
performance.v_sound_cc(i)=a_c_in;
performance.gamma_gc_t(i)=y.output.froz.gamma(2); % In throat
performance.Mm_gc_t(i)=Mm_gc_gola_in;    % In throat==chamber
performance.Tf_t(i)=Tf_gola_in;% In chamber
performance.rho_t(i)=rho_gola_in ;
performance.v_sound_t(i)=v_sound_t_in;
performance.v_t(i)=v_t_in;
performance.M_t(i)=Mach_t_in;

% 
% cea_performance.Isp(i)=Isp_cea;
% cea_performance.cf(i)=c_f_cea;
% cea_performance.cstar(i)=c_star_cea;

% Symbolic values to perform the calculation:
syms P_inc_Heox P_inc_Hef  T_inc_Heox T_inc_Hef V_inc_Heox V_inc_Hef ...
     V_inc_ox V_inc_f m_inc_ox m_inc_f P_inc_cc_ox P_inc_cc_f v_inc_ox ...
     v_inc_f y_inc x_inc M_esp_ox M_esp_f V_esp_ox V_esp_f mp_inc real

% Variables to do check:
Pc_min=20e5;
Mach_max=0.4;
RHe=R/Mm_He;

dt=0; % mettiamo questo per poi infittire solo a inizio e fine
while t(i)<=tmax && blowdown.Pc(i)>=Pc_min && blowdown.M_esp_tot_ox(i)<=Mi_ox && blowdown.M_esp_tot_f(i)<=Mi_fu % && blowdown.Pc(i)==blowdown.Pc_f && performance.M_cc(i)<=Mach_max 
    
    % Update the iteration index:
    i=i+1;
    dt = dt+i/2;
    if dt > 70
        dt = 50;
    end

    % if t(i-1) > 2400 % infittire alla fine
    %     dt = 10;
    % end
    

    % Update the time:
    t(i)=t(i-1)+dt;
 
    % The indeces are added after the solution of the system.
    % Write the system:
    eq1=M_esp_ox-dt*m_inc_ox;
    eq2=M_esp_f-dt*m_inc_f;
    
    eq3=V_esp_ox-M_esp_ox/rho_OX;
    eq4=V_esp_f-M_esp_f/rho_FU;
    
    eq5=V_inc_Heox-blowdown.V_Heox(i-1)-V_esp_ox;
    eq6=V_inc_Hef-blowdown.V_Hef(i-1)-V_esp_f;

    eq7=V_inc_ox-blowdown.V_ox(i-1)+V_esp_ox;
    eq8=V_inc_f-blowdown.V_f(i-1)+V_esp_f;
            
    eq9=blowdown.T_He_ox(i-1)*(blowdown.V_Heox(i-1))^(gamma_He-1)-T_inc_Heox*(V_inc_Heox)^(gamma_He-1);
    eq10=blowdown.T_He_f(i-1)*(blowdown.V_Hef(i-1))^(gamma_He-1)-T_inc_Hef*(V_inc_Hef)^(gamma_He-1);

    eq11=P_inc_Heox*V_inc_Heox-M_He_ox*RHe*T_inc_Heox;
    eq12=P_inc_Hef*V_inc_Hef-M_He_f*RHe*T_inc_Hef;

    eq13=m_inc_ox-cd_OX*A_ox_t*sqrt(2*rho_OX*x_inc*P_inc_cc_ox);
    eq14=m_inc_f-cD_fu*A_fu*sqrt(2*rho_FU*y_inc*P_inc_cc_f);

    eq15=v_inc_ox-m_inc_ox/(rho_OX*A_pipe_ox);
    eq16=v_inc_f-m_inc_f/(rho_FU*A_pipe_f);

    eq17=blowdown.P_Heox(i-1)-(1+x_inc)*P_inc_cc_ox-0.5*rho_OX*(v_inc_ox)^2-dP_feed;
    eq18=blowdown.P_Hef(i-1)-(1+y_inc)*P_inc_cc_f-0.5*rho_FU*(v_inc_f)^2-dP_feed;

    eq19=P_inc_cc_ox-P_inc_cc_f;
    
    eq20=mp_inc-m_inc_ox-m_inc_f; % conservation of mass 

    eq21=((mp_inc/(performance.rho_t(i-1)*At))/sqrt(performance.gamma_gc_t(i-1)*R/performance.Mm_gc_t(i-1)*performance.Tf_t(i-1)))-1;
    
    % Solve the system:
    sol=solve([eq1==0, eq2==0, eq3==0, eq4==0, eq5==0, eq6==0, eq7==0, ...
    eq8==0, eq9==0, eq10==0, eq11==0, eq12==0, eq13==0, eq14==0, ...
    eq15==0, eq16==0, eq17==0, eq18==0, eq19==0, eq20==0, eq21==0], ...
    [P_inc_Heox, P_inc_Hef,  T_inc_Heox, T_inc_Hef, V_inc_Heox, ...
    V_inc_Hef, V_inc_ox, V_inc_f, m_inc_ox, m_inc_f, P_inc_cc_ox, ...
    P_inc_cc_f,v_inc_ox, v_inc_f, x_inc, y_inc, M_esp_ox, M_esp_f, ...
    V_esp_ox, V_esp_f,mp_inc]);

    % Memorize the important variables:
    blowdown.P_Heox(i)=double(sol.P_inc_Heox);
    blowdown.P_Hef(i)=double(sol.P_inc_Hef);

    blowdown.T_He_ox(i)=double(sol.T_inc_Heox);
    blowdown.T_He_f(i)=double(sol.T_inc_Hef);

    blowdown.v_ox(i)=double(sol.v_inc_ox);
    blowdown.v_f(i)=double(sol.v_inc_f);

    blowdown.m_ox(i)=double(sol.m_inc_ox);
    blowdown.m_f(i)=double(sol.m_inc_f);
    blowdown.mp(i)=double(sol.mp_inc);

    blowdown.Pc(i)=double(sol.P_inc_cc_ox);
    blowdown.Pc_f(i)=double(sol.P_inc_cc_f);
   
    blowdown.x(i)=double(sol.x_inc);
    blowdown.y(i)=double(sol.y_inc);

    blowdown.V_Heox(i)=double(sol.V_inc_Heox);
    blowdown.V_Hef(i)=double(sol.V_inc_Hef);

    blowdown.V_ox(i)=double(sol.V_inc_ox);
    blowdown.V_f(i)=double(sol.V_inc_f);

    blowdown.M_esp_ox(i-1)=double(sol.M_esp_ox);
    blowdown.M_esp_f(i-1)=double(sol.M_esp_f);

    blowdown.V_esp_ox(i-1)=double(sol.V_esp_ox);
    blowdown.V_esp_f(i-1)=double(sol.V_esp_f);

    blowdown.M_esp_tot_ox(i)=blowdown.M_esp_tot_ox(i-1)+blowdown.M_esp_ox(i-1);
    blowdown.M_esp_tot_f(i)=blowdown.M_esp_tot_f(i-1)+blowdown.M_esp_f(i-1);

    % Post-processing:
    blowdown.r_of(i)=blowdown.m_ox(i)/blowdown.m_f(i);
    
    % vedere se possiamo inserire la dipendenza dalla temperatura in cea di
    % LOX e RP1
   %  y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',blowdown.r_of(i),'case','RP1_LOX','p,bar',blowdown.Pc(i)*1e-5,'sup',epsilon,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end');
     y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',blowdown.r_of(i),'case','RP1_LOX','fac','acat',eps_c,'p,bar',blowdown.Pc(i)*1e-5,'sup',epsilon,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end');

    % The output of CEA is selected as:
    %   1: injector
    %   2: comb.end
    %   3: Throat
    %   4: AR=200 
    % Select from CEA reasults:
    performance.Mm_gc(i)=y.output.froz.mw(1);       % In throat==chamber
    cp_in=(y.output.froz.cp_tran.froz(1)*1e3);
    performance.gamma_gc(i)=cp_in/(cp_in-(R/performance.Mm_gc(i)));
    performance.Tf(i)=y.output.froz.temperature(1); % In chamber
    performance.rho(i)=y.output.froz.density(1);    % In chamber

    performance.Mm_gc_t(i)=y.output.froz.mw(3);       % In throat==chamber
    cp_t=(y.output.froz.cp_tran.froz(3)*1e3);
    performance.gamma_gc_t(i)=cp_t/(cp_t-(R/performance.Mm_gc_t(i)));
    performance.Tf_t(i)=y.output.froz.temperature(3); % In chamber
    performance.rho_t(i)=y.output.froz.density(3);    % In chamber
    
    % cea performance are data used only to confront and not for further
    % computation
    % cea_performance.Isp(i)=y.output.froz.isp_vac(end);
    % cea_performance.cf(i)=y.output.froz.cf_vac(end);
    % cea_performance.cstar(i)=y.output.froz.cstar(end);
    
    % Calculate performances:
    f_Pe=@(Pe) (((performance.gamma_gc(i)+1)/2)^(1/(performance.gamma_gc(i)-1))) * ((Pe/blowdown.Pc(i))^(1/performance.gamma_gc(i))) *...
    sqrt(((performance.gamma_gc(i)+1)/(performance.gamma_gc(i)-1))*(1-(Pe/blowdown.Pc(i))^((performance.gamma_gc(i)-1)/performance.gamma_gc(i))))-(1/epsilon);
    Pe_guess=y.output.froz.pressure(4)*1e5;
    performance.Pe(i)=fsolve(f_Pe,Pe_guess,options);
    performance.ve(i)=sqrt(2*(performance.gamma_gc(i)/(performance.gamma_gc(i)-1))*(R/performance.Mm_gc(i))*performance.Tf(i)*(1-(performance.Pe(i)/blowdown.Pc(i))^((performance.gamma_gc(i)-1)/performance.gamma_gc(i))));
    performance.c_star(i)=sqrt((R/performance.Mm_gc(i)*performance.Tf(i))/(performance.gamma_gc(i)*(2/(performance.gamma_gc(i)+1))^((performance.gamma_gc(i)+1)/(performance.gamma_gc(i)-1))));

    performance.T(i)=blowdown.mp(i)*performance.ve(i)+Ae*performance.Pe(i);
    performance.Isp(i)=performance.T(i)/(9.81*blowdown.mp(i));
    performance.cf(i)= performance.T(i)/(blowdown.Pc(i)*At);
    performance.Itot=performance.Itot+performance.T(i)*dt;
    performance.T_2D(i)=lambda*blowdown.mp(i)*performance.ve(i)+Ae*performance.Pe(i);

    % velocità in camera
    performance.v_cc(i)=blowdown.mp(i)/(Ac*performance.rho(i)); % hyp: conservazione di massa
    performance.v_sound_cc(i)=sqrt(performance.gamma_gc(i)*R/performance.Mm_gc(i)*performance.Tf(i));
    performance.M_cc(i)=performance.v_cc(i)/performance.v_sound_cc(i);
    %velocità in gola
    performance.v_t(i)=blowdown.mp(i)/(At*performance.rho_t(i)); % hyp: conservazione di massa
    performance.v_sound_t(i)=sqrt(performance.gamma_gc_t(i)*R/performance.Mm_gc_t(i)*performance.Tf_t(i));
    performance.M_t(i)=performance.v_t(i)/performance.v_sound_t(i);

end

%% verify condition of exit
 if blowdown.M_esp_tot_ox(end)>=Mi_ox 
        disp('The mass of oxidizer expelled is bigger than the available')
    elseif blowdown.M_esp_tot_f(end)>=Mi_fu
        disp('The mass of fuel expelled is bigger than the available')
    elseif blowdown.Pc(end)~=blowdown.Pc_f
        disp('Boundary conditions of pressure on chamber pressure from fuel and oxidixer lines does not match')
    elseif performance.M_cc(end)>Mach_max 
      disp('The Mach number in combusiton chamber is not acceptable')
 end
 %%

eps1=(blowdown.Pc(1)*At/blowdown.mp(1))/(performance.c_star(1))
eps2=(blowdown.Pc(end)*At/blowdown.mp(end))/(performance.c_star(end))

% performances plot
figure

tiledlayout(3,2)


nexttile
hold on 
grid minor
plot(t(2:end),performance.T_2D(2:end),LineWidth=1.5)
ylim([0,1200])
xlabel('time-[s]')
ylabel('Thrust-[N]')
title('T(s)')

nexttile
hold on 
grid minor
plot(t(2:end),performance.Isp(1,2:end),LineWidth=1.5)
ylim([320 360])
xlabel('time-[s]')
ylabel('Isp-[s]')
title('Isp(s)')

nexttile
hold on 
grid minor
plot(t(2:end),blowdown.r_of(1,2:end),LineWidth=1.5)
ylim([2.1 2.9])
xlabel('time-[s]')
ylabel('O/F')
title('O/F(s)')

nexttile
hold on 
grid minor
plot(t(2:end),blowdown.Pc(1,2:end)*1e-5,LineWidth=1.5)
ylim([19 60])
xlabel('time-[s]')
ylabel('Pc-[bar]')
title('Pc(s)')

nexttile
hold on 
grid minor
plot(t(2:end),blowdown.mp(1,2:end),t(2:end),blowdown.m_ox(1,2:end),t(2:end),blowdown.m_f(1,2:end),LineWidth=1.5)
ylim([0 0.5])
xlabel('time-[s]')
ylabel('mp-[kg/s]')
legend('mp','mox','mf')
title('m(s)')

nexttile
hold on 
grid minor
plot(t(2:end),performance.c_star(1,2:end),LineWidth=1.5)
ylim([1500 2000])
xlabel('time-[s]')
ylabel('cstar-[m/s]')
title('cstar(s)')

figure

tiledlayout(3,2)


nexttile
hold on 
grid minor
plot(t(2:end),blowdown.P_Hef(2:end),LineWidth=1.5)
xlabel('time-[s]')
ylabel('Pressure-[N]')
title('P(s)')

nexttile
hold on 
grid minor
plot(t(2:end),blowdown.P_Heox(2:end),LineWidth=1.5)
xlabel('time-[s]')
ylabel('Pressure-[Pa]')
title('P(s)')

blowdown.M_esp_tot_ox(end) = capra(7);
nexttile
hold on 
grid minor
plot(t(2:end),blowdown.M_esp_tot_f(2:end),LineWidth=1.5)
yline(capra(8))
xlabel('time-[s]')
ylabel('Mass-[kg]')
title('M(s)')

nexttile
hold on 
grid minor
plot(t(2:end),blowdown.M_esp_tot_ox(2:end),LineWidth=1.5)
yline(capra(7))
xlabel('time-[s]')
ylabel('Mass-[kg]')
title('M(s)')


B_predicted=blowdown.P_Heox(1)/blowdown.P_Heox(end)
% %%
%     % Check conditions
%     if t(end) > tmax
%         disp('Condition: t(end) <= tmax is not met');
% 
%     elseif blowdown.Pc(end) < Pc_min
%         disp('Condition: blowdown.Pc(end) >= Pc_min is not met');
%     elseif blowdown.Mesp_ox(end) > Mi_ox
%         disp('Condition: blowdown.Mesp_ox(end) <= Mi_ox is not met');
%     elseif blowdown.Mesp_fu(end) > Mi_fu
%         disp('Condition: blowdown.Mesp_fu(end) <= Mi_fu is not met');
%     end
