function [Dv,r_dep,r_flyby,r_arr,Dv1_vect,Dv2_vect,Dv_pow,vt1_in,vt2_in,vt1_out,vt2_out,v_inf_minus,v_inf_plus] = dV_interplanetary_ToF(t1,ToF1,ToF2)
% DESCRIPTION:
%  Function to compute the total cost of the mission
% --------------------------------------------------------------------------
% PROTOTYPE:
%  [Dv,r_dep,r_flyby,r_arr,Dv1_vect,Dv2_vect,Dv_pow,vt1_in,vt2_in,vt1_out,
%  vt2_out,v_inf_minus,v_inf_plus] = dV_interplanetary_ToF(t1,ToF1,ToF2)
% --------------------------------------------------------------------------
% INPUT:
%  t1[1]          Departure date in MJD2000                                [days]
%  ToF1[1]        Time of flight of the first Lambert arc in MJD2000       [days]
%  ToF2[1]        Time of flight for the second Lambert arc in MJD2000     [days]
% --------------------------------------------------------------------------
% OUTPUT:
%  Dv[1]            DeltaV total cost of the mission                       [km/s]
%  r_dep[3x1]       State vector of the departure planet                   [km]
%  r_flyby[3x1]     State vector of the flyby planet                       [km]
%  r_arr[3x1]       State vector of the arrival object                     [km]
%  Dv1_vect[1]      DeltaV cost to enter in the First Lambert Arc          [km/s]
%  Dv2_vect[1]      DeltaV arrival cost                                    [km/s]
%  Dv_pow[1]        DeltaV cost of the powered gravity assist              [km/s]
%  vt1_in[3x1]      Initial velocity vector of the first Lambert arc 
%                   in Cartesian coordinates                               [km/s]
%  vt2_in[3x1]      Initial velocity vector of the second Lambert arc 
%                   in Cartesian coordinates                               [km/s]
%  vt1_out[3x1]     Final velocity vector of the first Lambert arc 
%                   in Cartesian coordinates                               [km/s]
%  vt2_out[3x1]     Final velocity vector of the second Lambert arc 
%                   in Cartesian coordinates                               [km/s]
%  v_inf_minus[3x1] Incoming excess velocity vector            
%                   in Cartesian coordinates                               [km/s]
%  v_inf_plus[3x1]  Outgoing excess velocity vector
%                   in Cartesian coordinates                               [km/s]
% --------------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% --------------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% --------------------------------------------------------------------------

%Sun gravitational constant
mu_Sun = astroConstants(4);

% Definition of the departure and arrival planets
dep_planet = 1; % Mercury
flyby_planet = 3; % Earth
arr_body = 21; % Asteroid N.21

Dv=NaN;
r_flyby=NaN;
r_arr=NaN;
Dv1_vect=NaN;
Dv2_vect=NaN;
Dv_pow=NaN;
vt1_in=NaN;
vt2_in=NaN;
vt1_out=NaN;
vt2_out=NaN;
v_inf_minus=NaN;
v_inf_plus=NaN;

if t1<date2mjd2000([2058, 01, 01, 00, 00, 00]) 
t2=t1+ToF1;
t3=t2+ToF2;
[kep_dep, ~]  = uplanet(t1, dep_planet);
[r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_Sun);
[kep_flyby, ~] = uplanet(t2, flyby_planet);

%First Lambert Arc:Mercury-Earth
[r_flyby,  v_flyby] = kep2car(kep_flyby(1), kep_flyby(2), kep_flyby(3), kep_flyby(4), kep_flyby(5), kep_flyby(6), mu_Sun);
[~, ~, ~, ERROR1, vt1_in, vt1_out, ~, ~] = lambertMR(r_dep, r_flyby, ToF1 * (24*3600), mu_Sun, 0, 0, 0, 2);
vt1_in = vt1_in(:); vt1_out = vt1_out(:);

%Second Lambert Arc: Earth-Asteroid N.21
kep_arr = ephNEO(t3, arr_body);
[r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_Sun);
[~,~,~,ERROR2, vt2_in, vt2_out,~,~] = lambertMR(r_flyby, r_arr,ToF2 * (24*3600), mu_Sun, 0, 0, 0, 2);
vt2_in = vt2_in(:); vt2_out = vt2_out(:);

if ERROR1 == 0 && ERROR2 == 0 
    v_inf_minus= vt1_out-v_flyby;
    v_inf_plus=vt2_in-v_flyby;
    if isfinite(v_inf_plus(1)) && isfinite(v_inf_minus(2))
    %Powered flyby cost
    [Dv_pow,r_p_hyp]=PoweredFlyBy(v_inf_minus,v_inf_plus);
        if not(isnan(r_p_hyp))
        %Cost to enter in the First Lambert Arc
        Dv1_vect=norm(vt1_in-v_dep);
        %Arrival cost
        Dv2_vect = norm(v_arr - vt2_out);
        Dvtot_vect = Dv1_vect + Dv2_vect+ Dv_pow;
        Dv = Dvtot_vect;
        end
    end
end
end
end

