function v_rot=rodrigues(v, u, del)
% DESCRIPTION:
% Function to rotate a vector v of an angle delta around unit vector u
% (counterclockwise)
% ------------------------------------------------------------------
% PROTOTYPE:
% v_rot=rodrigues(v, u, del)
% -------------------------------------------------------------------
% INPUT:
%   v[3x1]       vector to be rotated                        [km/s]
%   u[3x1]       vector around which the rotation happens    [-]
%   del[1]       rotation angle                              [rad]
% -------------------------------------------------------------------
% OUTPUT:
%   v_rot[3x1]   rotated vector                              [km/s]
% -------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% -------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% -------------------------------------------------------------------
v_rot=v*cos(del)+cross(u,v)*sin(del)+u*dot(u,v)*(1-cos(del));

