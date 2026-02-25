
function [ WMM_coeffs, K, g, h]=IGRF_coeffs(N)
% Compute the Schmidt quasi-normalization of IGRF coefficients up to 
%     order 13 and compute the matrix K useful to the on-line calculation of
%     the Gaussian normalized associated Legendre polynomials and their
%     derivatives.
% 
% 
% 
%     Notes: 
%     -All the matrices calculated are [n x m];
%     -g and h are in nT.
% 

%% Schmidt's quasi-normalization coefficients

S_nm = zeros(N,N+1); 
S_00 = 1;
% i=n j=m+1
for i = 1:N 
    if i==1
        S_nm(i,1) = S_00* (2*i-1)/i;
    else
        S_nm(i,1) = S_nm(i-1,1) * (2*i-1)/i;
    end
end

for i = 1:N
    for j = 1:N       
            S_nm(i,j+1) = S_nm(i,j)*sqrt( (i-j+1)*(kroneckerDelta(sym(j),1) + 1 )/ (i+j) );      
    end
end

%% g,h
WMM_coeffs = readmatrix("igrf13coeffs.txt");
WMM_coeffs = WMM_coeffs(4:end, [2 3 end-1]);

%N+1 because m need to go from 0 to N so N+1 elements
g = zeros(N,N+1);
h = zeros(N,N+1);

for i = 1:length(WMM_coeffs(:,1))
    if WMM_coeffs(i, 2) == 0
        g(WMM_coeffs(i,1),WMM_coeffs(i,2)+1) = WMM_coeffs(i,3);
    elseif WMM_coeffs(i, 2) == WMM_coeffs(i-1, 2)
        h(WMM_coeffs(i,1),WMM_coeffs(i,2)+1) = WMM_coeffs(i,3);
    else
    g(WMM_coeffs(i,1),WMM_coeffs(i,2)+1) = WMM_coeffs(i,3);
    end
end

%% Normalization
for i = 1:N
    for j = 1:N+1
        g(i,j) = S_nm(i,j)*g(i,j);
        h(i,j) = S_nm(i,j)*h(i,j);
    end
end

%% K

K = zeros(N,N+1);

for i = 2:N               
    for j = 1:i+1        
            K(i,j) = ((i-1)^2 - (j-1)^2) / ( (2*i-1) * (2*i-3) );
    end
end
