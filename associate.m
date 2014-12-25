% function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           Q                   2X2
%           z_i(t)              2X1
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1X1
%           outlier             1X1
%           nu^i(t)             2XN
%           S^i(t)              2X2XN
%           H^i(t)              2X3XN
function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)

N = size(M,2);
zhat = zeros(2,N);
H = zeros(2,3,N);
S = zeros(2,2,N);
nu = zeros(2,N);
D_m = zeros(N,1);
phi = zeros(N,1);

for j=1:N
    zhat(:,j) = observation_model(mu_bar, M,j);
    H(:,:,j) = jacobian_observation_model(mu_bar,M,j,zhat,j);
    S(:,:,j) = H(:,:,j)*sigma_bar*(H(:,:,j)') + Q;
    nu(:,j) = z_i - zhat(:,j);
    nu(2,j) = mod(nu(2,j)+pi,2*pi)-pi;
    D_m(j) = (nu(:,j)') * inv(S(:,:,j)) * (nu(:,j));
    phi(j) = det( 2*pi * S(:,:,j) )^(-1/2) * exp(-1/2 * D_m(j));
end

[m,c] = max(phi);
outlier = D_m(c)>=Lambda_m;

end