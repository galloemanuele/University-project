function [acc_pert] = acc_pert_fun_J2_moon(t, state, mu, J2, mu_M, R_E, frame)
%
% Function to compute the perturbing acceleration due to the J2 effect 
% and moon perturbations in cartesian coordinates and keplerian elements 
% (through Gauss's planetary equations)
% 
% PROTOTYPE:
%  [acc_pert]=acc_pert_fun_J2_moon(t, state, mu, J2, mu_M, R_E, frame)
% 
% INPUT:
%  t [1]          Time (MJD2000) expressed in seconds [s]
%  state [6,1]    Orbital state in cartesian coordinates state=[rx,ry,rz,vx,vy,vz] in XYZ case or keplerian in the others  [km],[km/s]
%  mu[1]          Gravitational parameter of the Earth   [km^3/s^2]
%  J2 [1]         Second zonal harmonic of the Earth   [-]
%  mu_M [1]       Gravitational parameter of the Moon   [km^3/s^2]
%  R_E [1]        Mean radius of the Earth    [km]
%  frame [str]    Choose the frame the perturbation must be evaluated between: 
%                 -'XYZ': ECI, Earth-Centred Inertial Reference Frame
%                 -'RSW': Radial--Transversal--Out-of-plane reference frame
%                 -'TNH': Tangential – Normal – Out-of-plane reference frame
%
% OUTPUT:
% acc_pert [3,1]  Perturbed acceleration due to J2 effect and moon pertubation [km/s^2]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version

if strcmp(frame, 'XYZ')
    r_XYZ=state(1:3); % position vector of the spacecraft [km]
    t_MJD2000=t/86400;              % time MJD2000 [days]
    [rm, ~]=ephMoon(t_MJD2000);     % extract the Moon position vector from the ephemerides [km]
    rmnorm=norm(rm);                % modulus of the Moon position vector [km]
    rms=rm'-r_XYZ;                  % position vector of the Moon wrt the s/c [km]
    rmsnorm=norm(rms);              % modulus of Moon position vector wrt the s/c [km]

    % Perturbation acceleration due to Moon gravity:
    a_moon = mu_M.*((rms./(rmsnorm^3))-(rm'./(rmnorm^3))); % [km/s^2]

elseif strcmp(frame, 'RSW')
    %Calculate the acceleration in XYZ frame:
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);
    [r_XYZ,v_XYZ]=kep2car(a,e,i,OM,w,f);
    t_MJD2000=t/86400;              % time MJD2000 [days]
    [rm, ~]=ephMoon(t_MJD2000);     % extract the Moon position vector from the ephemerides [km]
    rmnorm=norm(rm);                % modulus of the Moon position vector [km]
    rms=rm'-r_XYZ;                  % position vector of the Moon wrt the s/c [km]
    rmsnorm=norm(rms);              % modulus of Moon position vector wrt the s/c [km]

    % Perturbation acceleration due to Moon gravity:
    a_moon = mu_M.*((rms./(rmsnorm^3))-(rm'./(rmnorm^3))); % [km/s^2]
    %The perturbation of the moon in RSW reference frame is:
    a_moon = XYZ2RSW(a_moon, r_XYZ, v_XYZ, mu);


elseif strcmp(frame, 'TNH')
    %Calculate the acceleration in XYZ frame:
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);
    [r_XYZ,v_XYZ]=kep2car(a,e,i,OM,w,f);
    t_MJD2000=t/86400;              % time MJD2000 [days]
    [rm, ~]=ephMoon(t_MJD2000);     % extract the Moon position vector from the ephemerides [km]
    rmnorm=norm(rm);                % modulus of the Moon position vector [km]
    rms=rm'-r_XYZ;                  % position vector of the Moon wrt the s/c [km]
    rmsnorm=norm(rms);              % modulus of Moon position vector wrt the s/c [km]

    % Perturbation acceleration due to Moon gravity in XYZ frame:
    a_moon = mu_M.*((rms./(rmsnorm^3))-(rm'./(rmnorm^3))); % [km/s^2]
    %The perturbation of the moon in TNH reference frame is:
    a_moon = XYZ2TNH(a_moon, r_XYZ,v_XYZ);
else 
    error('The frame chosen is not correct')
end 

% Perturbing acceleration due to J2:
a_J2=acc_pert_fun_J2(t, state, mu, J2, R_E, frame); % [km/s^2]

% Total perturbing acceleration due to moon and J2:
acc_pert(1)=a_J2(1)+a_moon(1);
acc_pert(2)=a_J2(2)+a_moon(2);
acc_pert(3)=a_J2(3)+a_moon(3);

acc_pert=[acc_pert(1); acc_pert(2); acc_pert(3)];