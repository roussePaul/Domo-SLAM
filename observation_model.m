% function h = observation_model(x,j,psi)
% This function is the implementation of the h function.
% The bearing should lie in the interval [-pi,pi)
% Inputs:
%           x(t)        (3+2*N)X1
%           j           1X1
%           psi         1x1
% Outputs:  
%           h           1X1
function h = observation_model(x,j,psi)

[N,Ne,Nf,nf] = defSizes(x);

%% Parameters
d=0.2;


i=Nf+nf*(j-1);
rho = x(i);
t = x(i+1);
t = getAngle(x(4),t);
xr = x(1:(Nf-1));
h = -(xr(1)*cos(t)+xr(2)*sin(t) - rho + d*cos(xr(3)-t)) / cos(t - psi - xr(3));

end