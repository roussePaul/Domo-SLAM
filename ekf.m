% function [mu,sigma,outliers] = ekf(mu,sigma,R,Q,z,known_associations,u,M,Lambda_M,t)
% This function should perform one iteration of EKF localization.
% Inputs:
%           mu(t-1)             3X1
%           sigma(t-1)          3X3
%           R                   3X3
%           Q                   2X2
%           z                   2Xn
%           known_associations  1Xn
%           u                   3X1
%           M                   2XN
%           t                   1X1
%           Lambda_m            1X1
%           Map_IDS             1XN
% Outputs:
%           mu(t)               3X1
%           sigma(t)            3X3
%           outliers            1X1
function [mu,sigma,outliers] = ekf(mu,sigma,R,Q,z,known_associations,u,M,Lambda_m,Map_IDS,t)
[mu_bar,sigma_bar] = predict(mu,sigma,u,R);
n = size(z,2);
USE_KNOWN_ASSOCIATIONS = 1;
USE_BATCH_UPDATE = 0;
outliers = 0;
count = 0;

for i = 1 : n
    [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z(:,i),M,Lambda_m,Q);
    map_id = find(Map_IDS == known_associations(i));
    if c ~= map_id
        display(sprintf('warning, %d th measurement(of landmark %d) was incorrectly associated to landmark %d, t=%d',i,map_id,c,t));
    end
    if outlier
        display(sprintf('%d th measurement was labeled as outlier, t=%d',i,t));
        outliers = outliers + 1;
        continue;
    end
    if USE_KNOWN_ASSOCIATIONS
        c = map_id;
    end
    count = count + 1;
    nu_bar = squeeze(nu(:,c));
    S_bar = squeeze(S(:,:,c));
    H_bar = squeeze(H(:,:,c));
    [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar);   
end
mu = mu_bar;
sigma = sigma_bar;
    
if sum(sum(sigma~=sigma'))
    display('warning, sigma is not symmetric');
end
end