% function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           NeX1
%           sigma_bar(t)        NeXNe
%           Q                   1X1
%           z_i(t)              1X1
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1X1
%           outlier             1X1
%           nu^i(t)             1XN
%           S^i(t)              1X1XN
%           H^i(t)              1XNeXN
function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,angleMeasure,Lambda_m,Q)


Ne = size(mu_bar,1);

N = (Ne-3)/2;
zhat = zeros(N,1);
H = zeros(1,Ne,N);
S = zeros(N,1);
nu = zeros(N,1);
D_m = zeros(N,1);
phi = zeros(N,1);

for j=1:N
    zhat(j) = observation_model(mu_bar,j,angleMeasure);
    H(:,:,j) = jacobian_observation_model(mu_bar,j,zhat,j,angleMeasure);
    S(j) = H(:,:,j)*sigma_bar*(H(:,:,j)') + Q;
    nu(j) = z_i - zhat(j);
    D_m(j) = nu(j) * inv(S(j)) * nu(j);
    phi(j) = det( 2*pi * S(j) )^(-1/2) * exp(-1/2 * D_m(j));
end


[m,c] = max(phi);
outlier = D_m(c)>=Lambda_m;

end