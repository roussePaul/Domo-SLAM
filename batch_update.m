% function [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar)
% This function should perform the update process(single update).
% You need to make sure that the output sigma_bar is symmetric.
% The last line makes sure that ouput sigma_bar is always symmetric.
% Inputs:
%           mu_bar(t)       NeX1
%           sigma_bar(t)    NeXNe
%           H_bar(t)        nX3
%           Q_bar(t)		nXn
%           nu_bar(t)       2nX1
% Outputs:
%           mu(t)           3X1
%           sigma(t)        3X3
function [mu,sigma] = batch_update(mu_bar,sigma_bar,H_bar,Q_bar,nu_bar)
Ne = size(mu_bar,1);
K = sigma_bar * (H_bar') * inv( H_bar * sigma_bar * (H_bar') + Q_bar);
mu = mu_bar + K * nu_bar;
sigma = (eye(Ne)-K * H_bar) * sigma_bar;
sigma = (sigma + sigma')/2;
end