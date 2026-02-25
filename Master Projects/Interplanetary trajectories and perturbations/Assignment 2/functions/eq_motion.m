function [ds] = eq_motion(t, state, acc_pert, mu, frame)
% Function that evaluates the equations of motion (Cartesian or Keplerian).
%
% PROTOTYPE:
%  [ds] = eq_motion(t, state, acc_pert, mu, frame)
%
% INPUTS:
%  t [1]           Time (MJD2000) expressed in seconds [s]
%  state [6,1]     Orbital state in cartesian coordinates state=[rx,ry,rz,vx,vy,vz] in XYZ case or keplerian in the others  [km],[km/s]
%  acc_pert [3,1]  Perturbing acceleration considered [km/s^2]
%  mu[1]           Gravitational parameter of the Earth   [km^3/s^2]
%  frame [str]     Choose the frame the perturbation must be evaluated between: 
%                  -'XYZ': ECI, Earth-Centred Inertial Reference Frame
%                  -'RSW': Radial--Transversal--Out-of-plane reference frame
%                  -'TNH': Tangential – Normal – Out-of-plane reference frame
%
% OUTPUT:
%  ds [6,1]        derivative of orbital state in the frame chosen  [km/s],[km/s^2]
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
   r = state(1:3);
   v = state(4:6);
   rnorm = norm(r);
   ds = [v; -mu/(rnorm^3)*r + acc_pert];

% RSW: Radial - Transversal - Out-of-plane Reference Frame
elseif strcmp(frame, 'RSW')
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);

    ar = acc_pert(1);
    as = acc_pert(2);
    aw = acc_pert(3);

    % Auxiliary values
    b = a*sqrt(1-e^2);
    p = b^2/a;
    h = sqrt(p*mu);
    r = p/(1 + e*cos(f));
    vel = sqrt(2*mu/r - mu/a);

    % Derivatives
    da = 2*a^2/h*(e*sin(f)*ar + p*as/r);
    de = (p*sin(f)*ar + ((p+r)*cos(f) + r*e)*as)/h;
    di = r*cos(f+w)*aw/h;
    dOM = r*sin(f+w)*aw/(h*sin(i));
    dw = (-p*cos(f)*ar + (p+r)*sin(f)*as)/(h*e) - dOM*cos(i);
    df = h/(r^2) + (p*cos(f)*ar - (p+r)*sin(f)*as)/(e*h);

    ds = [da; de; di; dOM; dw; df];


% TNH: Tangential – Normal – Out-of-plane Reference Frame
elseif strcmp(frame,'TNH')
    a = state(1);
    e = state(2);
    i = state(3);
    OM = state(4);
    w = state(5);
    f = state(6);

    at = acc_pert(1);
    an = acc_pert(2);
    ah = acc_pert(3);

    % Auxiliary values
    b = a*sqrt(1-e^2);
    p = b^2/a;
    h = sqrt(p*mu);
    r = p/(1 + e*cos(f));
    vel = sqrt(2*mu/r - mu/a);

    % Derivatives
    da = 2*a^2*vel/mu*at;
    de = (2*(e + cos(f))*at - r*sin(f)*an/a)/vel;
    di = r*cos(f+w)*ah/h;
    dOM = r*sin(f+w)*ah/h/sin(i);
    dw = (2*sin(f)*at + (2*e + r*cos(f)/a)*an)/(e*vel) - dOM*cos(i);
    df = h/r^2 - (2*sin(f)*at + (2*e + r*cos(f)/a)*an)/(e*vel);

    ds = [da; de; di; dOM; dw; df];
end

end

