function [Dvi,Dvf,thf,Dt,at,et]=bitangent(ai,ei,wi,af,ef,wf,type,mu)
% DESCRIPTION:
%  Function to calculate the required parameters of the specified type of
%  bitangent transfer between the two orbits.
% -------------------------------------------------------------------------
% PROTOTYPE:
%  [Dvi,Dvf,thf,Dt,at,et]=bitangent(ai,ei,wi,af,ef,wf,type,mu)
% -------------------------------------------------------------------------
% INPUT:
%  ai[1]      Semi-major axis of initial orbit            [km]
%  ei[1]      Eccentricity of initial orbit               [-]
%  wi[1]      Periapsis argument of initial orbit         [rad]
%  af[1]      Semi-major axis of final orbit              [km]
%  ef[1]      Eccentricity of final orbit                 [-]
%  wf[1]      Periapsis argument of final orbit           [rad]
%  type[str]    Type of orbital transfer. Choose among 'pa', which stands for
%               apoapsis to periapsis transfer, 'ap',which stands for 
%               periapsis to apoapsis transfer, 'aa', which stands for 
%               apoapsis to apoapsis transfer, and 'pp',which stands for 
%               periapsis to periapsis transfer.
% mu[1]       Standard gravitational parameter            [km^3/s^2]
%               If not specified,mu default value is 
%               Earth standard gravitational parameter.
% -------------------------------------------------------------------------
% OUTPUT:
%  Dvi[1]     First impulse cost of the maneuver              [km/s]
%  Dvf[1]     Second impulse cost of the maneuver             [km/s]
%  thf[1]     Final true anomaly                              [rad]
%  Dt[1]      Time frame spent to accomplish the maneuver     [s]
%  at[1]      Semi-major axis of the transfer orbit           [km]
%  et[1]      Eccentricity of the transfer orbit              [-]
% -------------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% -------------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% -------------------------------------------------------------------------

if nargin==7
    mu=398600.433;
end

p_i=ai*(1-ei^2);
p_f=af*(1-ef^2);
rpi=p_i/(1+ei);
rai=p_i/(1-ei);
rpf=p_f/(1+ef);
raf=p_f/(1-ef);
err=0;
thf=0;
if mod(wf-wi,2*pi)==0 %eccentricity vectors are characterized by same verse
   if strcmp(type,'ap')==1
     ra_t_temp=rai; 
     rp_t_temp=rpf;
   elseif strcmp(type,'pa')==1
     ra_t_temp=raf;
     rp_t_temp=rpi;
   else 
       err=1;
       error('insert "ap" or "pa" otherwise change the eccentricities');
   end
elseif mod(abs(wf-wi),2*pi)==pi %eccentricity vectors are characterized by opposite verse
   if strcmp(type,'pp')==1
     rp_t_temp=rpi; 
     ra_t_temp=rpf;
   elseif strcmp(type,'aa')==1
     ra_t_temp=raf; 
     rp_t_temp=rai;
   else
       err=1;
       error('insert "aa" or "pp" otherwise change the eccentricities');
   end 
end 

if type(2)=='a'
    thf=pi;
end

%Check over ra_t and rp_t values. This check is justified by the 
%possibility to determine an orbit in which apoapsis and periapsis are 
%switched if compared to initial orbit parameters.
rp_t=min(ra_t_temp,rp_t_temp);
ra_t=max(ra_t_temp,rp_t_temp);
if rp_t~=rp_t_temp
    disp('Transfer orbit apoapsis and periapsis are switched in relation to initial ones');
end 

if err==0
    at=(ra_t+rp_t)/2;
    et=(ra_t-rp_t)/(ra_t+rp_t);
    Dvi=abs(sqrt(2*mu*(1/rp_t-1/(2*at)))-sqrt(2*mu*(1/rp_t-1/(2*ai))));
    Dvf=abs(sqrt(2*mu*(1/ra_t-1/(2*af)))-sqrt(2*mu*(1/ra_t-1/(2*at))));
    Dt=pi*sqrt((at^3)/mu);
else
    [Dvi,Dvf,thf,Dt]=zeros(1,4);
end