function [T,Y]=orbitPropagator(y0,tspan,choice,mu)
% DESCRIPTION:
% Function to propagate the orbit
% ----------------------------------------------------------------------
% PROTOTYPE:
% [T,Y]=orbitPropagator(y0,tspan,choice,mu)
% ----------------------------------------------------------------------
% INPUT:
%  y0[6x1]     State of the body (rx, ry, rz, vx, vy, vz)        [km, km/s]          
%  tspan[1xn]  Vector of times                                   [s]
%  choice:     'up' for unperturbed 2 body problem
%              'p' for perturbed 2 body problem
%  mu[1]       Planetary constant of the planet (mu = mass * G)  [km^3/s^2]
% ----------------------------------------------------------------------
% OUTPUT:
%  T[nx1]      Time vector of the state                          [s]
%  Y[nx6]      State of the body ( rx, ry, rz, vx, vy, vz )      [km, km/s]
% ----------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% ----------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% ----------------------------------------------------------------------
if nargin==3
   mu = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
end
if strcmp(choice,'up')
        % Set options for the ODE solver
        options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
        % Perform the integration
        [ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
elseif strcmp(choice,'p')
        par = astroConstants([23,9])
        R_e=par(1);
        J2=par(2);
        % Set options for the ODE solver
        options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
        % Perform the integration
        [ T, Y ] = ode113( @(t,y) ode_2bp_p(t,y,mu,R_e,J2), tspan, y0, options );
else 
    error("insert p to choose perturbed 2bp and up to choose unperturbed one")
end

end