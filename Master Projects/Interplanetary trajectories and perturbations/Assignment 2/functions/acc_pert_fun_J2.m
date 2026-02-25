function [acc_pert] = acc_pert_fun_J2(t, state, mu, J2, R, frame)
%
% Function to compute the perturbing acceleration due to the J2 effect in 
% cartesian coordinates and keplerian elements (through Gauss's planetary 
% equations).
% 
% PROTOTYPE:
%  [acc_pert] = acc_pert_fun_J2(t, state, mu, J2, R, frame)
% 
% INPUT:
%  t [1]          Time (MJD2000) expressed in seconds [s]
%  state [6,1]    Orbital state in cartesian coordinates state=[rx,ry,rz,vx,vy,vz] in XYZ case or keplerian in the others  [km],[km/s]
%  mu[1]          Gravitational parameter of the Earth   [km^3/s^2]
%  J2 [1]         Second zonal harmonic of the Earth   [-]
%  R [1]          Mean radius of the Earth    [km]
%  frame [str]    Choose the frame the perturbation must be evaluated between: 
%                 -'XYZ': ECI, Earth-Centred Inertial Reference Frame
%                 -'RSW': Radial--Transversal--Out-of-plane reference frame
%                 -'TNH': Tangential – Normal – Out-of-plane reference frame
%
% OUTPUT:
% acc_pert [3,1]  Perturbed acceleration due to J2 effect [km/s^2]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
% 09/01/2024: First version

if strcmp(frame, 'XYZ')
    r = state(1:3);
    v = state(4:6);
    rnorm = norm(r);
    x = r(1);
    y = r(2);
    z = r(3);
    % Acceleration due to J2 in Cartesian position
    acc_pert = 3/2*((J2*mu*R^2)/rnorm^4)*[(x/rnorm)*(5*(z^2/rnorm^2)-1);
                                          (y/rnorm)*(5*(z^2/rnorm^2)-1);
                                          (z/rnorm)*(5*(z^2/rnorm^2)-3)];

    % RSW: Radial - Transversal - Out-of-plane Reference Frame
elseif strcmp(frame, 'RSW')
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);

    [r, ~] = kep2car(a, e, i, OM, w, f, mu);
    rnorm = norm(r);

    % Acceleration due to J2 in Radial - Transversal - Out-of-plane position
    acc_pert = -3/2*((J2 *mu*R^2)/rnorm^4).*[1-3*(sin(i)^2)*(sin(f+w))^2; 
                                            (sin(i)^2)*sin(2*(f+w)); 
                                            sin(2*i)*sin(f+w)];
    
    % TNH: Tangential – Normal – Out-of-plane Reference Frame
    elseif strcmp(frame,'TNH')
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);

    % Auxiliary values
    b = a*sqrt(1-e^2);
    p = b^2/a;
    h = sqrt(p*mu);
    r = p/(1 + e*cos(f));
    vel = sqrt(2*mu/r - mu/a);
    
    [r, ~] = kep2car(a, e, i, OM, w, f, mu);
    rnorm = norm(r);

    % Acceleration due to J2 in Radial - Transversal - Out-of-plane position
    acc_pert = -3/2*((J2 *mu*R^2)/rnorm^4).*[1-3*(sin(i)^2)*(sin(f+w))^2; 
                                            (sin(i)^2)*sin(2*(f+w)); 
                                            sin(2*i)*sin(f+w)];
    ar = acc_pert(1);
    as = acc_pert(2);
    aw = acc_pert(3);

    % Acceleration due to J2 in Tangential – Normal – Out-of-plane position
    % Rotation from Radial - Transversal - Out-of-plane position
    A = h/(p*vel)*[e*sin(f), -(1+e*cos(f)); 1+e*cos(f), e*sin(f)];
    b = [ar; as];
    x = A\b;
    at = x(1);
    an = x(2);
    ah = aw;
        
    acc_pert = [at; an; ah];

end