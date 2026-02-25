function [r,v] = kep2car(a, e, i, OM, om, th, mu)
% Conversion from Keplerian elements to Cartesian coordinates.
%
% PROTOTYPE:
%  [r,v] = kep2car(a, e, i, OM, om, th, mu)
%
% INPUT:
%  a [1x1]          Semi-major axis [km]
%  e [1x1]          Eccentricity [-]
%  i [1x1]          Inclination [rad]
%  OM [1x1]         RAAN [rad]
%  om [1x1]         Pericentre anomaly [rad]
%  th [1x1]         True anomaly [rad]
%  mu [1x1]         Standard gravitational parameter [km^3/s^2]. 
%
% OUTPUT:
%  r [3x1]          Position vector [km]
%  v [3x1]          Velocity vector [km/s]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version

if nargin == 6
   mu=398600.433;
end

%Calculate p and r:
p=a*(1-e^2);
r=p/(1+e*cos(th));

%Find the status vectors in the perifocal system:
r_pf=r*[cos(th);sin(th);0];
v_pf=sqrt(mu/p)*[-sin(th);e+cos(th);0];

%Define the rotation arrays:
R3=@(a)[cos(a) sin(a) 0;
    -sin(a) cos(a) 0;
    0 0 1];
R1=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

%Use the rotation arrays to determine the status vectors in cartesian
%system:
r=R3(OM)'*R1'*R3(om)'*r_pf;
v=R3(OM)'*R1'*R3(om)'*v_pf;