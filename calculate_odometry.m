% function u = calculate_odometry(e_l,e_R,E_T,B,delta_t,mu)
% This function should calculate the odometry information
% Inputs:
%           e_L(t):         1X1
%           e_R(t):         1X1
%           E_T:            1X1
%           B:              1X1
%           R_L:            1X1
%           R_R:            1X1
%           delta_t:        1X1
%           mu(t-1):        3X1
% Outputs:
%           u(t):           3X1
function u = calculate_odometry(delta_t,mu,v,w)
if ~delta_t
    u = [0;0;0];
    return;
end
u = delta_t * [v*cos(mu(3));v*sin(mu(3));w];
end