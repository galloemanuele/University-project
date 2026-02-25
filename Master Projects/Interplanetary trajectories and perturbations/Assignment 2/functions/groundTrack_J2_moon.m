function [alpha,delta,lon,lat]=groundTrack_J2_moon(type, par, thetaG0, t_vect, mu_E, J2, mu_M, R_E, w_E)
% Function that calculates the right ascension, the declination, the
% longitude and latitude, useful for the ground track in the perturbed case
% due to J2 effect and moon perturbation.
% 
% PROTOTYPE:
%  [alpha,delta,lon,lat]=groundTrack_J2_moon(type, par, thetaG0, t_vect, mu_E, J2, mu_M, R_E, w_E)
% 
% INPUT:
%  type [char]               Describes how the orbit initial state is given ("kep" = Keplerian Elements; "car" = Cartesian Coordinates)
%  par [6,1]                 Orbit initial state vector which with Keplerian Elements [km],[rad] in case of Keplerian or initial position and velocity in case of 'Cartesian' [km],[km/s]
%  thetaG0 [1]               Right ascension of Greenwich meridian at t0    [rad]
%  t_vect [1,N]              Time span vector (N is the number of points used for the linspace) [s]
%  mu_E [1]                  Gravitational parameter of the Earth [km^3/s^2]
%  J2 [1]                    Second zonal harmonic [-]
%  mu_M [1]                  Gravitational parameter of the Moon   [km^3/s^2]
%  R_E [1]                   Mean radius of the Earth    [km]
%  w_E [1]                   Earth's spin rate [rad/s]. If not assigned, it's automatically set to 7.291597763887421e-05 [rad/s]
%
% OUTPUT: 
%  alpha_deg [N,1]                       vector of right-ascension value at each time-step [rad]
%  delta_deg [N,1]                       vector of declination value at each time-step [rad]
%  lon_wrapped [N,1]                     vector of longitude value at each time-step [deg]
%  lat [N,1]                             vector of latitude value at each time-step [deg]
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
    
if strcmp(type,"kep")
    [R0,V0] = kep2car(par(1),par(2), par(3), par(4), par(5), par(6), mu_E);
elseif  strcmp(type,"car")
    R0=par(1:3);
    V0=par(4:6);
else
    error("insert 'car' or 'par' in the field type")
end

y0=[R0,V0];
% Set options for the ODE solver:
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration:
acc_pert = @(t,s) acc_pert_fun_J2_moon(t, s, mu_E, J2, mu_M, R_E, 'XYZ'); % perturbed acceleration
ds = @(t,s) eq_motion(t, s, acc_pert(t,s), mu_E, 'XYZ'); % resulting perturbing state
[T, Y] = ode113(ds, t_vect, y0, options);

% Matrix containing R and V
R=Y(:,1:3);

% Preallocate vectors of unknown for speed:
delta=zeros(1,length(T));
alpha=zeros(1,length(T));

for i=1:length(T)
    normR=norm(R(i,:));

    % Declination:
    delta(i)=asin(R(i,3)/normR);

    % Right ascension:
    if R(i,2)/normR>0
        alpha(i)=acos(R(i,1)/(normR*cos(delta(i))));
    else
        alpha(i)=2*pi-acos(R(i,1)/(normR*cos(delta(i))));
    end
end

%True anomaly of Greenwich meridian at a generic time t:
t_0=ones(1,length(T))*t_vect(1);
thetaG=thetaG0*ones(1,length(T))+w_E.*(T'-t_0); 

%Longitude:
lon=wrapTo([-pi,pi],alpha-thetaG);
lon=rad2deg(lon);

%Latitude:
lat=wrapTo([-pi/2,pi/2],delta); 
lat=rad2deg(lat);

end