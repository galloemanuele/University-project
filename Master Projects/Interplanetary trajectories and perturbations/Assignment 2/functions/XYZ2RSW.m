function x_RSW = XYZ2RSW(x_XYZ, r_XYZ, v_XYZ, mu)
% function to convert from Cartesian coordinates to tangential-out of
% plane-normal reference frame
% 
% PROTOTYPE:
%  x_RSW = XYZ2RSW(x_XYZ, r_XYZ,v_XYZ)
%  
% INPUT:
%  x_XYZ [3,1]    Generic quantity whose components are in inertial frame [undefined a priori]
%  r_XYZ [3,1]    Position vector in inertial frame [km]
%  v_XYZ [3,1]    Velocity vector in inertial frame [km/s]
% 
% OUTPUT:
%  x_RSW [3,1]    Generic quantity whose components are in RSW frame [undefined a priori]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version

% Rotate the vector x_XYZ in TNH frame:
x_TNH=XYZ2TNH(x_XYZ, r_XYZ,v_XYZ);

% Conversion in keplerian elements:
[a, e, ~, ~, ~, f] = car2kep(r_XYZ, v_XYZ, mu);

% Calculate useful quantities:
b = a*sqrt(1-e^2);
p = b^2/a;
h = sqrt(p*mu);
r = p/(1 + e*cos(f));
vel = sqrt(2*mu/r - mu/a);

% Rotation matrix from TNH to RSW:
A=h/(p*vel)*[e*sin(f), -(1+e*cos(f)); 1+e*cos(f), e*sin(f)];

% Rotate the vector x_TNH in RSW frame:
x_RS=A*[x_TNH(1); x_TNH(2)];
x_W=x_TNH(3);

x_RSW=[x_RS; x_W];

