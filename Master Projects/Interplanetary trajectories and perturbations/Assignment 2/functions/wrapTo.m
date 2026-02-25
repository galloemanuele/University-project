function [out]=wrapTo(interval,in)
% Function to wrap the angles specified in input in the interval chosen.
%
% PROTOTYPE:
%  [out]=wrapTo(interval,in)
%
% INPUT:
%  in [N.D., 1]     input vector of the angles that must be wrapped. "N.D"
%                   stands for undefined, since it can be used an arbitrary dimension
%                   vector for the input. [rad]
%  interval [1,2]   interval in which the function wraps the input angles.
%                   the first element must be negative. [rad, rad]
% OUTPUT:
%  out [N.D., 1]    vector of angles wrapped. "N.D" stands for undefined, 
%                   since it can be used an arbitrary dimension vector for
%                   the output. [rad]
%
% CONTRIBUTORS:
%  Casiero Alessia
%  Fiume Elisa 
%  Gallo Emanuele 
%  Marotta Arianna
% 
% VERSIONS:
%  09/01/2024: First version
for i=1:length(in)
    
    while in(i)>interval(2)
        in(i)=in(i)-abs(interval(2)+abs(interval(1)));
    end
    while abs(in(i))>abs(interval(1)) && in(i)<0
        in(i)=in(i)+abs(interval(2)+abs(interval(1)));
    end
end
out=in;
end