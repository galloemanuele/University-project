function P=BuildGrid(t2,dep_min,arr_max,T_dep,T_arr)
% Function that builds the grid of local minima point for n flyby dates
% PROTOTYPE:
% P=BuildGrid(t2,dep_min,arr_max,T_dep,T_arr)
% --------------------------------------------------------------------
% INPUT:
%  t2[nx1]          Flyby chosen date vector in MJD2000      [Days]
%  dep_min[1]       Minimum departure date  in MJD2000       [Days]
%  arr_max[1]       Maximum departure date  in MJD2000       [Days]
%  T_dep[1]         Period of the departure planet           [Days]
%  T_arr[1]         Period of the arrival object             [Days]
% OUTPUT:
% P[x,y,n]          3D matrix contatining the local minima   [Days]
%                   coordinates(dates) 
%                   x and y cannot be defined a priori
% --------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% --------------------------------------------------------------------
% Fixed flyby date vector considering only the chosen time window for
% departure
t1=linspace(dep_min,dep_min+T_dep,30);
for i=1:length(t2)
    t3=linspace(t2(i),t2(i)+T_arr,30);
    for j=1:length(t1)
    for k=1:length(t3)
        Dvtot_mat(j,k,i) = dV_interplanetary_unc(t1(j),t2(i),t3(k));
    end
    end
end
%Calculate the minimum for each rectangle 
t=zeros(length(t2),3);
Dv_min=zeros(length(t2),1);
for i=1:length(t2)
    t3=linspace(t2(i),t2(i)+T_arr,100);
    [Dv_min(i),k_opt]=min(min(Dvtot_mat(:,:,i)));
    [~,j_opt]=min(Dvtot_mat(:,k_opt,i));
    t(i,:)=[t1(j_opt) t2(i) t3(k_opt)];
    
end

%Creation of the grid of estimated minima point position
x=zeros(length(t2),length(t(end,1):T_dep:t(end,2)));
y=zeros(length(t2),length(t(end,2):T_arr:arr_max));

for i=1:length(t2)
    n1=length(t(i,1):T_dep:t(i,2));
    n2=length(t(i,3):T_arr:arr_max);
    x(i,1:n1)=t(i,1):T_dep:t(i,2);
    y(i,1:n2)=t(i,3):T_arr:arr_max;
end
x(x==0)=NaN;
y(y==0)=NaN;

%Make the grid for each flyby date
for i=1:length(t2)
   pos=0;
   for j=1:length(x(i,:))
       for k=1:length(y(i,:))
           P(k+pos,:,i)=[x(i,j) y(i,k)];
       end
       pos=pos+k;
   end
end