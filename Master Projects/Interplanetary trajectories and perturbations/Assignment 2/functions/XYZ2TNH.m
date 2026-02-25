function x_TNH = XYZ2TNH(x_XYZ, r_XYZ,v_XYZ)
% function to convert from Cartesian coordinates to tangential-out of
% plane-normal reference frame
% 
% PROTOTYPE:
%  x_TNH = XYZ2TNH(x_XYZ, r_XYZ,v_XYZ)
%  
% INPUT:
%  x_XYZ [3,1]    Generic quantity whose components are in inertial frame [undefined a priori]
%  r_XYZ [3,1]    Position vector in inertial frame [km]
%  v_XYZ [3,1]    Velocity vector in inertial frame [km/s]
% 
% OUTPUT:
%  x_TNH [3,1]    Generic quantity whose components are in TNH frame [undefined a priori]
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version

% Define the unit vectors composing the TNH frame:
t = v_XYZ/norm(v_XYZ);                             % Tangential versor
h = cross(r_XYZ,v_XYZ)/norm(cross(r_XYZ,v_XYZ));   % Out of plane versor
n = cross(h,t);                                    % Normal versor

% The rotation matrix is:
A = [t,n,h]';

% The quantity rotated from inertial to TNH frame is:
x_TNH=A*x_XYZ;