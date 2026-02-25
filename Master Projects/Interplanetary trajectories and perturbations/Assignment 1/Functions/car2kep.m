function kep = car2kep(r, v, mu)
%PROTOTYPE:
% kep = car2kep(r, v, mu)
%--------------------------------------------------------------------------
%DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements.
%--------------------------------------------------------------------------
%INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Standard gravitational parameter [km^3/s^2]. If not specified,
%          mu default value is earth standard gravitational parameter.
%--------------------------------------------------------------------------
% OUTPUT:
% kep[6x1] Keplerian elements vector [km,-,rad,rad,rad,rad]

%Calculate the absolute value of status vectors:
r_norm=norm(r);
v_norm=norm(v);

%Specific angular moment and its absolute value:
h=cross(r,v);
h_norm=norm(h);

%Inclination of the orbit:
i=acos(h(3,1)/h_norm);

%Eccentricity vector and eccentricity:
e_vett=1/mu*((v_norm^2-mu/r_norm)*r-(dot(r,v))*v);
e=norm(e_vett);

%Specific mechanical energy:
E=1/2*v_norm^2-mu/r_norm;
%Semi-major axe:
a=-mu/(2*E);

%Nodes line:
K=[0 0 1]';
N=cross(K,h);
N_norm=norm(N);

%RAAN of the orbit:
if(N(2,1)>=0)
    OM=acos(N(1,1)/N_norm);
else
    OM=2*pi-acos(N(1,1)/N_norm);
end

%Periapsis argument:
if(e_vett(3,1)>=0)
    om=acos(dot(N,e_vett)/(N_norm*e));
else
    om=2*pi-acos(dot(N,e_vett)/(N_norm*e));
end

%Radial velocity:
vr=dot(r,v)/r_norm;

%True anomaly:
if(vr>=0)
    th=acos(dot(e_vett,r)/(e*r_norm));
else
    th=2*pi-acos(dot(e_vett,r)/(e*r_norm));
end
kep=[a, e, i, OM, om, th];
end