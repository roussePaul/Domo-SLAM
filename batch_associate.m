% function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           (3+2*Ne)X1
%           sigma_bar(t)        (3+2*Ne)X(3+2*Ne)
%           Q                   1X1
%           z(t)                1Xn
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1Xn
%           outlier             1Xn
%           nu_bar(t)           nX1
%           H_bar(t)            nX3
function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,angle,Lambda_m,Q)

[N,Ne,Nf,nf,sM] = defSizes(mu_bar);

n = size(z,2);

zhat = zeros(N,n);
H = zeros(sM,Ne,N);
S = zeros(N,1);

nu = zeros(N,1);
D_m = zeros(N,1);
phi = zeros(N,1);
c = zeros(n,1);
outlier = zeros(n,1);


nu_bar=[];
H_bar=[];

for i=1:n
    for j=1:N
        zhat(j) = observation_model(mu_bar,j,angle(i));
        H(:,:,j) = jacobian_observation_model(mu_bar,j,zhat,j,angle(i));
        S(j) = H(:,:,j)*sigma_bar*(H(:,:,j)') + Q;
        
        nu(j) = z(i) - zhat(j);
        D_m(j) = (nu(j)') * inv(S(j)) * nu(j);
        phi(j) = det( 2*pi * S(j) )^(-1/2) * exp(-1/2 * D_m(j));
    end
    [m,c(i)] = max(phi);
    outlier(i) = D_m(c(i))>=Lambda_m;
    
    nu_bar = [nu_bar nu(c(i))'];
    H_bar = [H_bar H(:,:,c(i))'];
    
    %nu_f(:,i) = nu(:,c(i));
    %H_f(:,:,i) = H(:,:,c(i))';
end

%index = find(outlier==0);
%np = size(index,1);
%H_bar = H_f(:,:,index);
%H_bar = reshape(H_bar,3,np*2);
%nu_bar = nu_f(:,index);

%nu_bar = reshape(nu_bar,1,2*np);

nu_bar = nu_bar';
H_bar = H_bar';



end