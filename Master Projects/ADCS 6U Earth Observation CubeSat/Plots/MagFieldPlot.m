clc 
clear 
close
R_E=6371;
environment.R = R_E; 
settings.N=13;
[WMM, environment.WMM.K, environment.WMM.g, environment.WMM.h]=IGRF_coeffs(settings.N);
N=13;
w_E=7.291597763887421e-05;

long=linspace(-pi,pi,1000)
lat=linspace(-pi/2,pi/2,1000)
r=500+R_E;
alpha_G=0;
B=zeros(length(long),length(lat))
for i=1:length(long)
    for j=1:length(lat)
    
    B(i,j) = norm(magneticField(alpha_G,long(i),r,lat(j),environment,N));

    end
end
B_val=max(norm(B))

%%
[C,h]=contour(rad2deg(long),rad2deg(lat),B',20000:2000:52000,'LineWidth', 2,'ShowText','on','HandleVisibility','off');
clabel(C,h,'FontSize',5,'Color','black')
hcb = colorbar;
hcb.Limits = [2*10^4, 5.2*10^4];
set(get(hcb,'Title'),'String','B [nT]');
hold on
I=imread('EarthTexture_half_trans.jpg');
X=[-180,180];
Y=[90,-90];
h=image(X,Y,I);
uistack(h,'bottom')
xlabel('Longitude $[^\circ]$', 'Interpreter', 'Latex')
ylabel('Latitude $[^\circ]$', 'Interpreter', 'Latex')

%%
function B = magneticField(alpha_G,long,radius,lat,environment,N)
%{  
    Return the magnetic field component in inertial reference frame

%}
theta = (pi/2-lat);     % theta is the colatitude [rad]
alpha = long + alpha_G; 

g = environment.WMM.g;
h = environment.WMM.h;
K = environment.WMM.K;
R = environment.R;


%% P-dP
% validated through the P and dP tables in "Jeremy Davis, Mathematical Modeling of Earth's Magnetic
% Field, 2004"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_{(n),(m+1)}: 
% given P(a,b) -> n=a, m=b-1
% n from 1 to N
% m from 0 to N
% N order of accuracy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(N,N+1);
dP = zeros(N,N+1);

P_00 = 1;
dP_00 = 0;

P_01 = 0; 
dP_01 = 0;

P(1,1) = cos(theta); % n=1, m=0
dP(1,1) = -sin(theta)*P_00; 

P(1,2) = sin(theta);  % n=1, m=1
dP(1,2) = cos(theta)*P_00;

for i = 2:N % ciclo in n a partire da n=2 a n=N
    for j = 0:i % ciclo in m da m=0 a m=n
        if j == i
        P(i,j+1) = sin(theta)*P(i-1,j);
        dP(i,j+1) = sin(theta)*dP(i-1,j) + cos(theta)*P(i-1,j);
        elseif i < 3
            if j == 0 % n=2, m=0
            P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P_00;
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1) ...
                - K(i,j+1)*dP_00;
            elseif j == 1 % n=2, m=1
                P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P_01;
                dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1)...
                    - K(i,j+1)*dP_01;
            end
        else
            P(i,j+1) = cos(theta)*P(i-1,j+1) - K(i,j+1)*P(i-2,j+1);
            dP(i,j+1) = -sin(theta)*P(i-1,j+1) + cos(theta)*dP(i-1,j+1) ...
                - K(i,j+1)*dP(i-2,j+1);
        end
    end
end


B_r = 0;
B_th = 0;
B_phi = 0;

%indice i scorre n da 1 a N
%indice j scorre m da 0 a n

for m = 0:N
    for n = 1:N
        if m<=n
            B_r = B_r + (R/radius)^(n+2)*(n+1)*...
                ((g(n,m+1)*cos(m*long) + h(n,m+1)*sin(m*long))*P(n, m+1));
            B_th = B_th + (R/radius)^(n+2)*...
                ((g(n,m+1)*cos(m*long) + h(n,m+1)*sin(m*long))*dP(n, m+1));
            B_phi = B_phi + (R/radius)^(n+2)*...
                (m*(-g(n,m+1)*sin(m*long) + h(n,m+1)*cos(m*long))* P(n, m+1));

        end

    end
end

B_th = -B_th;
B_phi = -B_phi/sin(theta);

B_x = ( B_r*cos(lat) + B_th*sin(lat) ) * cos(alpha) - B_phi*sin(alpha);
B_y = ( B_r*cos(lat) + B_th*sin(lat) ) * sin(alpha) + B_phi*cos(alpha);
B_z = B_r*sin(lat) - B_th*cos(lat);

B = [B_x B_y B_z]';
end