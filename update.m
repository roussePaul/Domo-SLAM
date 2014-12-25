% function [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar)
% This function should perform the update process(sequential update).
% You need to make sure that the output sigma_bar is symmetric.
% The last line makes sure that ouput sigma_bar is always symmetric.
% Inputs:
%           mu_bar(t)       NeX1
%           sigma_bar(t)    NeXNe
%           H_bar(t)        1XNe
%           S_bar(t)        1X1
%           nu_bar(t)       1X1
% Outputs:
%           mu_bar(t)       NeX1
%           sigma_bar(t)    NeXNe
function [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar)
Ne = size(mu_bar,1);
K = sigma_bar * (H_bar') * inv(S_bar);
mu_bar = mu_bar + K * nu_bar;
sigma_bar = (eye(Ne)- K*H_bar) * sigma_bar;
sigma_bar = (sigma_bar + sigma_bar')/2;
end