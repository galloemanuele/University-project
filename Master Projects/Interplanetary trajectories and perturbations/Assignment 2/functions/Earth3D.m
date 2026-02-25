function Earth3D(Rt)
% This function plots the Earth in the center of the 3D plot with the
% radius Rt specified
%
% PROTOTYPE:
%  Earth3D(Rt)
% 
% INPUT:
%  Rt [1]     Radius to be assigned to the Earth plot [km]
%             If not specified, Rt is automatically set to 6378.1363 km. 
%
% OUTPUT:
%  Plot
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
% 09/01/2024: First version

if nargin==0
    Rt=6378.1363;
end

load('topo.mat','topo','topomap1');
contour(0:359,-89:90,topo,[0 0],'b')
axis equal
hold on
image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
colormap(topomap1)
cla reset
axis equal
[x,y,z] = sphere(100);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
surface(x*Rt,y*Rt,z*Rt,props);
light('position',[1 1 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);