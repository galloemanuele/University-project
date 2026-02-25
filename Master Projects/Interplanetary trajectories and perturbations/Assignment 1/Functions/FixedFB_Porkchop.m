function FixedFB_Porkchop(t2,dep_min,arr_max,T_arr,T_dep)
% DESCRIPTION:
% Function to plot the graph showing the effectivness of the local minima
% grid method localization by comparing it with an actual porkchop plot with
% the same flyby date.
% ------------------------------------------------------------------
% PROTOTYPE:
% FixedFB_Porkchop(t2,dep_min,arr_max,T_arr,T_dep)
% -------------------------------------------------------------------
% INPUT:
%  t2[1]            Flyby chosen date in MJD2000    [Days]
%  dep_min[1]       Minimum departure date          [Days]
%  arr_max[1]       Maximum departure date          [Days]
%  T_arr[1]         Period of the arrival object    [Days]
%  T_dep[1]         Period of the departure planet  [Days]
% -------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% -------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% 

%% Porkchop plot with fixed flyby date
t1=linspace(dep_min,t2,200);
t3=linspace(t2,arr_max,200);
Dvtot_mat1=zeros(length(t1),length(t2));
for i=1:length(t1)
    i
    for j=1:length(t3)
        Dvtot_mat1(i,j) = dV_interplanetary_unc(t1(i),t2,t3(j));
    end
end

[Dv_min,j_opt]=min(min(Dvtot_mat1));
[~,i_opt]=min(Dvtot_mat1(:,j_opt));
t1_vect=zeros(1,length(t1));
t3_vect=zeros(1,length(t3));
for i=1:length(t1)
    t1_vect(i)=datenum(mjd20002date(t1(i)));
end
for i=1:length(t3)
    t3_vect(i)=datenum(mjd20002date(t3(i)));
end

h = figure(1);
hold on
contour(t1_vect,t3_vect,Dvtot_mat1' ,Dv_min+(0:1:30),'LineWidth', 1,'ShowText','off','HandleVisibility','off');
scatter3(t1_vect(i_opt),t3_vect(j_opt),Dv_min,'filled','o','MarkerFaceColor','r')
xlabel('Departure date')
ylabel('Arrival date')
datetick('x','yyyy mmm dd','keeplimits')
xlim([min(t1_vect), max(t1_vect)]);
datetick('y','yyyy mmm dd','keeplimits')
ylim([min(t3_vect), max(t3_vect)]);
grid on
hold on
hcb = colorbar;
hcb.Limits = [Dv_min, round(Dv_min)+30];
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
title("Original Porkchop plot with fixed FlyBy date")


h = figure(2);
hold on
contour(t1_vect,t3_vect,Dvtot_mat1' ,Dv_min+(0:1:30),'LineWidth', 1,'ShowText','off','HandleVisibility','off');
scatter3(t1_vect(i_opt),t3_vect(j_opt),Dv_min,'filled','o','MarkerFaceColor','r')
xlabel('Departure date')
ylabel('Arrival date')
datetick('x','yyyy mmm dd','keeplimits')
xlim([min(t1_vect), max(t1_vect)]);
datetick('y','yyyy mmm dd','keeplimits')
ylim([min(t3_vect), max(t3_vect)]);
grid on
hold on
hcb = colorbar;
hcb.Limits = [Dv_min, round(Dv_min)+30];
set(get(hcb,'Title'),'String','\DeltaV [km/s]');
title("Modified Porkchop plot with minima position estimation")
hold on


%% Initial minimum estimation: test to compare with the Porkchop plot
% Fixed flyby date

%Definition of arrival and departure dates restricted to the first
%"rectangle" of the pattern repetition

t1=linspace(dep_min,dep_min+T_dep,100);

for i=1:length(t2)
    t3=linspace(t2(i),t2(i)+T_arr,100);
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
    i
    t3=linspace(t2(i),t2(i)+T_arr,100);
    [Dv_min(i),k_opt]=min(min(Dvtot_mat(:,:,i)));
    [~,j_opt]=min(Dvtot_mat(:,k_opt,i));
    t(i,:)=[t1(j_opt) t2(i) t3(k_opt)];
end

%Minimum grid creation
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
P=[];
pos=0;
   for j=1:length(x(i,:))
       for k=1:length(y(i,:))
           
           P(k+pos,:)=[x(i,j) y(i,k)];
           if isnan(y(i,k)) || isnan(x(i,j))
               P_graph(k+pos,:)=[NaN NaN];
           else
               P_graph(k+pos,:)=[datenum(mjd20002date(x(i,j)))  datenum(mjd20002date(y(i,k)))];
           end
       end
       pos=pos+k;
   end
   
%Plot of the grid
for i=1:length(P_graph)
    hold on
    scatter3(P_graph(i,1),P_graph(i,2),dV_interplanetary_unc(P(i,1),t2,P(i,2)),'filled','o','MarkerFaceColor','r')
end

end
