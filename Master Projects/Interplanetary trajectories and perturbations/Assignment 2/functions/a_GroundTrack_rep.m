function [a]=a_GroundTrack_rep(type, m, k, mu, e, i, J2, R, w_E)
% 
% Function that computes the semi-major axis required to obtain a repeating
% ground track with the given values of k,m both in unperturbed and
% in the perturbed case, where the secular effects of J2 are taken into
% considerations.
% 
% PROTOTYPE:
%  [a]=a_GroundTrack_rep(type, m, k, mu, e, i, J2, R, w_E)
%
% INPUT:
% type [str]              Choose between 'perturbed' and 'unperturbed' case
%  k [1]                  Number of revolutions of the s/c to obtain the ground track repetition [-]
%  m [1]                  Number of revolutions of the planet to obtain the ground track repetition [-]
%  mu [1]                 Gravitational parameter of primary body [km^3/s^2]
%  e [1]                  S/c orbit's eccentricity [-]
%  i [1]                  S/c orbit's inclination [rad]
%  J2 [1]                 Second zonal harmonic [-]
%  R [1]                  Mean radius of the primary body [km]
%  w_E [1]                Spin rate of the primary body [rad/s]
%
% OUTPUT: 
%  a [1]                  Semi-major axis for the ground track repetition for the chosen case [km]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version

if nargin==8
    w_E=7.291597763887421e-05; %[rad/s]
end 

a_unp=((mu*m^2)/(w_E*k)^2)^(1/3);

if strcmp(type, 'unperturbed') % unperturbed case
    a=a_unp;
elseif strcmp(type, 'perturbed') % perturbed case
    OMdot = @(x) -(3/2*sqrt(mu)*J2*R^2/(1-e^2)^2/x^(7/2))*cos(i);
    omdot = @(x) -(3/2*sqrt(mu)*J2*R^2/(1-e^2)^2/x^(7/2))*(5/2*(sin(i))^2-2);
    M0dot = @(x) -(3/2*sqrt(mu)*J2*R^2/(1-e^2)^2/x^(7/2))*(1-3/2*(sin(i)^2));

    n = @(x) sqrt(mu/x^3);

    f = @(x) m/k-(w_E-OMdot(x))/(n(x) + omdot(x) + M0dot(x));
   
    options = optimset('TolX',1e-14);
    a = fzero(f, a_unp, options);
else 
    error('Select only "unperturbed" or "perturbed" type')
end 