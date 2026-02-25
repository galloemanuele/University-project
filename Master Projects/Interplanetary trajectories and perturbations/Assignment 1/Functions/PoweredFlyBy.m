function [Dv_pg,r_p,a_hyp1,e_hyp1,a_hyp2,e_hyp2,v_p1,v_p2]=PoweredFlyBy(v_inf_1,v_inf_2)
% DESCRIPTION:
% Function to compute the powered gravity assist manoeuvre around Earth.
% -------------------------------------------------------------------------
% PROTOTYPE:
% [Dv_pg,r_p,a_hyp1,e_hyp1,a_hyp2,e_hyp2,v_p1,v_p2]=PoweredFlyBy(v_inf_1,v_inf_2)
% ------------------------------------------------------------------------- 
% INPUT:
%  v_inf_1 [3,1]    incoming excess velocity vector                    [km/s]
%  v_inf_2 [3,1]    outgoing excess velocity vector                    [km/s]
% -------------------------------------------------------------------------
% OUTPUT:
%  Dv_pg [1]        delta v cost of the powered gravity assist         [km/s]
%  r_p [1]          radius of perigee of the flyby trajectory          [km]
%  a_hyp1[1]        incoming hyperbola semi-major axis                 [km]
%  a_hyp2 [1]       outgoing hyperbola semi-major axis                 [km]
%  e_hyp1 [1]       incoming hyperbola eccentricity                    [-]
%  e_hyp2 [1]       outgoing hyperbola eccentricity                    [-]
%  v_p1 [1]         incoming velocity magnitude at perigee             [km/s]
%  v_p2 [1]         outgoing velocity magnitude at perigee             [km/s]
%  v_inf_m_n [1]    magnitude of incoming excess velocity              [km/s]
%  v_inf_p_n [1]    magnitude of outgoing excess velocity              [km/s]
% -------------------------------------------------------------------------
% CONTRIBUTORS:
% Casiero Alessia
% Fiume Elisa
% Gallo Emanuele
% Marotta Arianna
% -------------------------------------------------------------------------
% VERSIONS:
%  01-01-2023: First version
% -------------------------------------------------------------------------

mu_E=astroConstants(13);
R_E=astroConstants(23);
%Height of Earth's atmosphere
h_atm=200; %[km]
%Turning angle
del=acos(dot(v_inf_2,v_inf_1)/(norm(v_inf_1)*norm(v_inf_2)));

%Solve the nonlinear system to find rp
rp_fun=@(r_p)  asin(1./(1+(r_p*norm(v_inf_1)^2)/mu_E))+asin(1./(1+(r_p*norm(v_inf_2)^2)/mu_E))-del;
options = optimset('Display','off','TolX',1e-14);
r_p=fzero(rp_fun,R_E+2000,options);

%Check if the radius obtained is physically feasible
r_p_limit=R_E+h_atm; %[km]
    if r_p>r_p_limit

    %Solution of the 2D hyperbola 1(in)
    v_inf_1_norm=norm(v_inf_1);
    a_hyp1=-mu_E/v_inf_1_norm^2;
    e_hyp1=(1+ r_p * v_inf_1_norm^2/mu_E);
    v_p1=mu_E/(sqrt(a_hyp1*mu_E*(1-e_hyp1^2)))*(1+e_hyp1);
    
    %Solution of the 2D hyperbola 2(out)
    v_inf_2_norm=norm(v_inf_2);
    a_hyp2=-mu_E/v_inf_2_norm^2;
    e_hyp2=(1+ r_p * v_inf_2_norm^2/mu_E);
    v_p2=mu_E/(sqrt(a_hyp2*mu_E*(1-e_hyp2^2)))*(1+e_hyp2);
    
    %Delta_v_p and Delta_v_flyby
    Dv_pg=abs(v_p1-v_p2);
    else
        Dv_pg=NaN;
        r_p=NaN;
    end
end
