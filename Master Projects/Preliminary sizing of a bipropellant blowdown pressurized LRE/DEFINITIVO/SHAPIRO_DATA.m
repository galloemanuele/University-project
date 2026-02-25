function [Re,Reff,C_f,Ldiv_id,gamma_eff]=SHAPIRO_DATA(eps_max)
R=8314.5; % J/K kmol
frozen=1; % at end of Combustion Chamber
eps_c = 10;
Pc=50;
o_f=2.30;
T_in_real=1e3; % data
losses_noz=1-0.9925; % hp
T_in_ideal=T_in_real*(1+losses_noz);
options = optimoptions('fsolve', 'Display', 'off');
eps=linspace(1,eps_max,100);

for i=1:length(eps)
 %y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',o_f,'case','RP1_LOX','p,bar',Pc_in,'sup',epsilon,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end');

 y=CEA('problem','rocket','frozen','nfz',frozen,'o/f',o_f,'case','RP1_LOX','fac','acat',eps_c,'p,bar',Pc,'sup',eps(i),'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'t(k)',300,'h,cal/mol',-5430,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end');

 % gamma_cea=y.output.froz.gamma(2)
% gamma_in=y.output.froz.gamma(2)*-y.output.froz.dlvpt(1) % gamma cea corretto con comprimibilità
 Mm_in=y.output.froz.mw(2); % throat
cp_in=(y.output.froz.cp_tran.froz(1)*1e3);
gamma_maggi=cp_in/(cp_in-(R/Mm_in));
gamma_in = gamma_maggi;

Tf_in=y.output.froz.temperature(1);

Pc_in=50e5; % bar, data
f=@(Pe) (((gamma_in+1)/2)^(1/(gamma_in-1))) * ((Pe/Pc_in)^(1/gamma_in)) *...
    sqrt(((gamma_in+1)/(gamma_in-1))*(1-(Pe/Pc_in)^((gamma_in-1)/gamma_in)))-(1/eps(i));

 Pe_guess=y.output.froz.pressure(4)*1e5;
 Pe_in=fsolve(f,Pe_guess,options);
 ve_in=sqrt(2*(gamma_in/(gamma_in-1))*(R/Mm_in)*Tf_in*(1-(Pe_in/Pc_in)^((gamma_in-1)/gamma_in)));
 c_star_in=sqrt((R/Mm_in*Tf_in)/(gamma_in*(2/(gamma_in+1))^((gamma_in+1)/(gamma_in-1))));
 
 At=T_in_ideal/((Pc_in*ve_in/c_star_in)+eps(i)*Pe_in);
 % mp_in=Pc_in*At/c_star_in;
 Ae=eps(i)*At;
 De=sqrt(4*Ae/pi);
 Reff(i)=De/2;
 Dt=sqrt(4*At/pi);

 alpha = deg2rad(15);
 Ldiv_id(i)=0.5*(De-Dt)/tan(alpha);
 Re(i)=y.output.froz.density(4)*ve_in*De/(y.output.froz.viscosity(4)*1e-6);
 C_f(i)=16/Re(i);
  Mm_in=y.output.froz.mw(4); % throat
 cp_in=(y.output.froz.cp_tran.froz(4)*1e3);
 gamma_eff(i)=cp_in/(cp_in-(R/Mm_in));
end
end