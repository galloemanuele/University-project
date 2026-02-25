% performances plot
clc
 clear all
% close all
%%
load("risultati_veri_al_3_sigma.mat")

t_nom=t;
T_nom=performance.T_2D(2:end);
Isp_nom=performance.Isp(1,2:end);
r_of_nom=blowdown.r_of(1,2:end);
Pc_nom=blowdown.Pc(1,2:end)*1e-5;
mp_nom=blowdown.mp(1,2:end);
m_ox_nom=blowdown.m_ox(1,2:end);
m_f_nom=blowdown.m_f(1,2:end);
c_star_nom=performance.c_star(1,2:end);
mach_nom=performance.M_t(1,2:end);

%%
load("arch_new_additive_minus.mat")

t_min=t;
T_min=performance.T_2D(2:end);
Isp_min=performance.Isp(1,2:end);
r_of_min=blowdown.r_of(1,2:end);
Pc_minimo=blowdown.Pc(1,2:end)*1e-5;
mp_min=blowdown.mp(1,2:end);
m_ox_min=blowdown.m_ox(1,2:end);
m_f_min=blowdown.m_f(1,2:end);
c_star_min=performance.c_star(1,2:end);
mach_min=performance.M_t(1,2:end);

%%
load("arch_new_additive_plus.mat")

t_plus=t;
T_plus=performance.T_2D(2:end);
Isp_plus=performance.Isp(1,2:end);
r_of_plus=blowdown.r_of(1,2:end);
Pc_plus=blowdown.Pc(1,2:end)*1e-5;
mp_plus=blowdown.mp(1,2:end);
m_ox_plus=blowdown.m_ox(1,2:end);
m_f_plus=blowdown.m_f(1,2:end);
c_star_plus=performance.c_star(1,2:end);
mach_plus=performance.M_t(1,2:end);


%%
figure

tiledlayout(3,2)


nexttile
hold on 
grid minor
plot(t_nom(2:end),T_nom,t_min(2:end),T_min,t_plus(2:end),T_plus,LineWidth=1.5)
ylim([0,1200])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('Thrust-[N]')
title('T(s)')

nexttile
hold on 
grid minor
plot(t_nom(2:end),Isp_nom,t_min(2:end),Isp_min,t_plus(2:end),Isp_plus,LineWidth=1.5)
ylim([340 360])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('Isp-[s]')
title('Isp(s)')

nexttile
hold on 
grid minor
plot(t_nom(2:end),r_of_nom,t_min(2:end),r_of_min,t_plus(2:end),r_of_plus,LineWidth=1.5)
ylim([1.5 3])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('O/F')
title('O/F(s)')

nexttile
hold on 
grid minor
plot(t_nom(2:end),Pc_nom,t_min(2:end),Pc_minimo,t_plus(2:end),Pc_plus,LineWidth=1.5)
ylim([19 60])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('Pc-[bar]')
title('Pc(s)')

nexttile
hold on 
grid minor
plot(t_nom(2:end),mp_nom,t_min(2:end),mp_min,t_plus(2:end),mp_plus,LineWidth=1.5)
ylim([0 0.5])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('mp-[kg/s]')
title('m(s)')

nexttile
hold on 
grid minor
plot(t_nom(2:end),c_star_nom,t_min(2:end),c_star_min,t_plus(2:end),c_star_plus,LineWidth=1.5)
ylim([1700 1800])
legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('cstar-[m/s]')
title('cstar(s)')
%%
figure
hold on 
grid minor
plot(t_nom(2:end),mach_nom,t_min(2:end),mach_min,t_plus(2:end),mach_plus,LineWidth=1.5)
ylim([0.8 1.2])
 legend('nom','minus','plus','Location','best')
xlabel('time-[s]')
ylabel('Mach')
title('Mach')

