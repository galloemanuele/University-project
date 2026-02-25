% iteration zero con incertezze additive -60 micron
clear all
close all
clc

format compact
uncertainty=+60e-6;
%%
T_in_real=1e3; % data
losses_noz=1-0.9925; % hp
T_in_ideal=T_in_real*(1+losses_noz);

epsilon=200; % hp

Pc_in=50; % bar, data
o_f=2.30; 
R=8314.5; % J/K kmol
frozen=1; % at end of Combustion Chamber
eps_c = 10;

 %y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',o_f,'case','RP1_LOX','p,bar',Pc_in,'sup',epsilon,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end');

 y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',o_f,'case','RP1_LOX','fac','acat',eps_c,'p,bar',Pc_in,'sup',epsilon,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end','screen');


% gamma_cea=y.output.froz.gamma(2)
% gamma_in=y.output.froz.gamma(2)*-y.output.froz.dlvpt(1) % gamma cea corretto con comprimibilità
Mm_in=y.output.froz.mw(2); % throat
cp_in=(y.output.froz.cp_tran.froz(1)*1e3);
gamma_maggi=cp_in/(cp_in-(R/Mm_in));
gamma_in = gamma_maggi;

Tf_in=y.output.froz.temperature(1);
Isp_cea=y.output.froz.isp_vac(end);
c_f_cea=y.output.froz.cf_vac(end);
c_star_cea=y.output.froz.cstar(end);
rho_in=y.output.froz.density(1);

gamma_gc_gola_in=y.output.froz.gamma(3); % In throat
Mm_gc_gola_in=y.output.froz.mw(3);       % In throat==chamber
Tf_gola_in=y.output.froz.temperature(3); % In chamber
rho_gola_in=y.output.froz.density(3);    % In chamber

%%
% simulation 
Pc_in=50e5; % bar, data
f=@(Pe) (((gamma_in+1)/2)^(1/(gamma_in-1))) * ((Pe/Pc_in)^(1/gamma_in)) *...
    sqrt(((gamma_in+1)/(gamma_in-1))*(1-(Pe/Pc_in)^((gamma_in-1)/gamma_in)))-(1/epsilon);

 Pe_guess=y.output.froz.pressure(4)*1e5;
 Pe_in=fsolve(f,Pe_guess);
 ve_in=sqrt(2*(gamma_in/(gamma_in-1))*(R/Mm_in)*Tf_in*(1-(Pe_in/Pc_in)^((gamma_in-1)/gamma_in)));
 c_star_in=sqrt((R/Mm_in*Tf_in)/(gamma_in*(2/(gamma_in+1))^((gamma_in+1)/(gamma_in-1))));
 
 At=T_in_ideal/((Pc_in*ve_in/c_star_in)+epsilon*Pe_in);
 mp_in=Pc_in*At/c_star_in;
 Ae=epsilon*At;
 De=sqrt(4*Ae/pi);
 Dt=sqrt(4*At/pi);

 I_sp=T_in_ideal/(9.81*mp_in);
 c_f=T_in_ideal/(Pc_in*At);

%% dim CC
L_star = 1.27;  % worst case scenario (?)% tbd range from 5 to 15
Ac = eps_c*At;
Dc = sqrt(4*Ac/pi);
Vc = L_star*At;
Lc = Vc/Ac;

% Mc verification
rho_c = Pc_in / (R/Mm_in * Tf_in); 
u_c_in = mp_in/(rho_in*Ac);
a_c_in =  sqrt(gamma_in*R/Mm_in*Tf_in);
M_cc_in =  u_c_in/a_c_in;




% Mc=0.2; % hypotesis
% Ac=At/Mc*((2/(gamma_in+1))*(1+((gamma_in-1)/2)*Mc^2))^((gamma_in+1)/(2*(gamma_in-1)));
% Dc=sqrt(4*Ac/pi);
% L_star=1.02; % hypotesis
% Vc=At*L_star;
% Lc=L_star*At/Ac;
% eps_c=Ac/At;
%% design of converegent nozzle and conical divergent nozzle

beta=deg2rad(45);
Lconv=0.5*(Dc-Dt)/tan(beta);

alpha = deg2rad(15);
Ldiv_id=0.5*(De-Dt)/tan(alpha);


%% maximum nozzle length

theta_n=deg2rad(28);
theta_e=deg2rad(4.9);

lambda=0.5*(1+cos((alpha+theta_e)/2));

T=lambda*mp_in*ve_in+Pe_in*Ae;
% [X,R]=bell_shape_madda(Dt,epsilon,28,4.9, Lc, Dc,1);

% alpha 2D in questo caso è 0.0075

%% minimum nozzle length

% Ldiv_min=0.6*Ldiv_id;
% theta_n=deg2rad(37);
% theta_e=deg2rad(13);
% 
% alpha_1=atan2(De-Dt,2*Ldiv_min);
% lambda=0.5*(1+cos((alpha_1+theta_e)/2));
% 
% T=lambda*mp_in*ve_in+Pe_in*Ae;

% [X,Raooo]=bell_shape_madda(Dt,epsilon,37,13, Lc, Dc,0.6);

%alpha 2D in questo caso è 0.0261

%% injection plate
clc
m_OX = (o_f/(1 + o_f))*mp_in;
m_FU = mp_in - m_OX;
rho_OX = 1450;
rho_FU = 870;
d_p = 0.25*Pc_in;
cD_fu = 0.7;
d_fu = 5e-4;
A_fu = m_FU/(cD_fu*sqrt(2*d_p*rho_FU));
A_inj_fu = pi*d_fu^2/4;
N_fu = A_fu/A_inj_fu
N_fu=ceil(N_fu);
A_fu_new = A_inj_fu*N_fu;
m_FU_new = A_fu_new*cD_fu*sqrt(2*d_p*rho_FU)

%% triplet
N_ox_t = 2*N_fu;
d_p_f = 0.15 * 2e6;
d_ox_t = 0.47e-3;
A_inj_ox_t = pi*d_ox_t^2/4;
A_ox_t = N_ox_t * A_inj_ox_t;
cd_OX = m_OX/(A_ox_t*sqrt(2*d_p*rho_OX))

m_fu_f_t = A_fu_new*sqrt(2*rho_FU*d_p_f);
m_ox_f_t = A_ox_t*sqrt(2*rho_OX*d_p_f);
m_ox_f_t/m_fu_f_t
m_ox_f_t+m_fu_f_t
% m_OX/m_FU_new
%  m_OX+m_FU_new

u_ox = cd_OX*sqrt(d_p/rho_OX);
u_fu = cD_fu*sqrt(d_p/rho_FU);
u_ox = cd_OX*sqrt(d_p_f/rho_OX);
u_fu = cD_fu*sqrt(d_p_f/rho_FU);

% control L_imp
L_D = 6;
L_imp = d_fu * L_D;
%% UPDATE VALUE OF DIAMETERS OF INJECTORS DUE TO ADDITIVE MANUFACTURING
% CASE 1: UNCERTAINTY OF DIAMETER
% CASE 1.1: - 60 micron
% uncertainty=-60e-6;
% CASE 1.2: 60 micron
% plus_uncertainty=60e-6;

d_ox_t = 0.47e-3+uncertainty;
A_inj_ox_t = pi*d_ox_t^2/4;
A_ox_t = N_ox_t * A_inj_ox_t;

d_fu = 5e-4+uncertainty;
A_fu = m_FU/(cD_fu*sqrt(2*d_p*rho_FU));
A_inj_fu = pi*d_fu^2/4;


%% check area pipes and velocity
clc
%  Area pipes
v0=15; % m/s

A_pipe_ox=m_OX/(rho_OX*v0);
D_pipe_ox=sqrt(4*A_pipe_ox/pi); % m

A_pipe_f=m_FU/(rho_FU*v0);
D_pipe_f=sqrt(4*A_pipe_f/pi);

% final velocity
of_f=2.2814;
m_f=0.2661;

m_OX_f = (of_f/(1 + of_f))*m_f;
m_FU_f = m_f - m_OX_f;

v_ox_f=m_OX_f/(rho_OX*A_pipe_ox);
v_fuel_f=m_FU_f/(rho_FU*A_pipe_ox);


%% pressurisation
dP_feed=50e3; % Pa
dP_inj_in=0.25*Pc_in;

P_tank_ox_in=Pc_in+0.5*rho_OX*v0^2+dP_feed+dP_inj_in;
P_inj_ox_in=Pc_in+0.5*rho_OX*v0^2+dP_feed;

P_tank_f_in=Pc_in+0.5*rho_FU*v0^2+dP_feed+dP_inj_in;
P_inj_f_in=Pc_in+0.5*rho_FU*v0^2+dP_feed;

%% maximum blowdown ratio
Pc_f=20e5; % Pa
dP_inj_f=0.25*Pc_f;

P_tank_ox_f=Pc_f+0.5*rho_OX*v_ox_f^2+dP_feed+dP_inj_in;
P_tank_f_f=Pc_f+0.5*rho_FU*v_fuel_f^2+dP_feed+dP_inj_in;

B_ox_max=P_tank_ox_in/P_tank_ox_f
B_f_max=P_tank_f_in/P_tank_f_f

%% estimation of tank volume
clc

B_ox_max=2.1;
B_f_max=2.1;
Ti_He_ox=90; % K
Ti_He_fu=330; % K
Dmax=1;
Hmax=2;
Vtot=Hmax*pi*(Dmax/2)^2;
Vtanks=Vtot*0.8-Vc;
gamma_He=1.67;

% ratio is the volumetric ratio between the tanks
 ratio=[1.36:0.001:1.39];

syms Vi_He_ox  Vi_He_fu Vf_He_ox  Vf_He_fu  Tf_He_ox Tf_He_fu  real
% eq 1, eq 2 Oxigen
% eq 3, eq 4 RP-1


for i=1:length(ratio)

    eq1=B_ox_max*Vi_He_ox^gamma_He-Vf_He_ox^gamma_He;
    eq2=Ti_He_ox*Vi_He_ox^(gamma_He-1)-Tf_He_ox*Vf_He_ox^(gamma_He-1);
    
    eq3=B_f_max*Vi_He_fu^gamma_He-Vf_He_fu^gamma_He;
    eq4=Ti_He_fu*Vi_He_fu^(gamma_He-1)-Tf_He_fu*Vf_He_fu^(gamma_He-1);
    
    eq5=Vf_He_ox-Vf_He_fu*ratio(i);
    eq6=Vf_He_ox+Vf_He_fu-Vtanks;
    
    % Risolvi il sistema di equazioni
    sol = solve([eq1==0,eq2==0,eq3==0,eq4==0,eq5==0,eq6==0],...
        [Vi_He_ox  Vi_He_fu Vf_He_ox  Vf_He_fu  Tf_He_ox Tf_He_fu]);
    
    Vi_ox=sol.Vf_He_ox-sol.Vi_He_ox;
    Mi_ox=Vi_ox*rho_OX;

    Vi_fu=sol.Vf_He_fu-sol.Vi_He_fu;
    Mi_fu=Vi_fu*rho_FU;

    r_of_test=Mi_ox/Mi_fu;

    GHIRINGHELLI(i,:)=[sol.Vi_He_ox  sol.Vi_He_fu sol.Vf_He_ox  sol.Vf_He_fu  sol.Tf_He_ox...
        sol.Tf_He_fu Mi_ox Mi_fu r_of_test];

end

GHIRINGHELLI=double(GHIRINGHELLI);

j=1;

for i=1:length(ratio)
    if(GHIRINGHELLI(i,end)<2.4 && GHIRINGHELLI(i,end)>2.2 )
        brandonisio(j,:)=[GHIRINGHELLI(i,:) ratio(i)];
        j=j+1;
    end
end
brandonisio

% brandonisio è l array con i risultati che rispettano il vincolo su o/f
% capra è la riga di brandonisio per cui o/f=2.29, valor medio tra massimo
% e minimo
%%
clc
 capra=[0.4649    0.3409    0.7250    0.5315   66.8296  245.0419  377.0887  165.8748    2.2733    1.3640];


% results from iteration zero
Vtank_ox=capra(3);
Vtank_fu=capra(4);

Tf_He_ox_prevista=capra(5);
Tf_He_fu_prevista=capra(6);


costante=capra(10);% rapoorto tra Vtank_ox/Vtank_f

% valori iniziali iterazione zero
Vi_He_ox=capra(1);
Vi_ox=Vtank_ox-Vi_He_ox;
Vi_He_fu=capra(2);
Vi_fu=Vtank_fu-Vi_He_fu;
Mi_ox=capra(7);
Mi_ox/(m_OX*60)
Mi_fu=capra(8);
Mi_fu/(m_FU*60)
ri_of_test=capra(9);
Mm_He=4.002602; % kg/kmol

M_He_ox=P_tank_ox_in*Vi_He_ox/(Ti_He_ox*R/Mm_He);
M_He_f=P_tank_f_in*Vi_He_fu/(Ti_He_fu*R/Mm_He);
M_He=M_He_f+M_He_ox;
v_cc_in=mp_in/(Ac*rho_in);
v_t_in=mp_in/(At*rho_gola_in);
v_sound_t_in=sqrt(gamma_gc_gola_in*R/Mm_gc_gola_in*Tf_gola_in);
Mach_t_in=v_t_in/v_sound_t_in;


save("RISULTATI_ITERAZIONE0.mat")

