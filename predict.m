% function [mu_bar,sigma_bar] = predict(mu,sigma,u,R)
% This function should perform the prediction step.
% Inputs:
%           mu(t-1)           3X1   
%           sigma(t-1)        3X3
%           u(t)              3X1
%           R                 3X3
% Outputs:   
%           mu_bar(t)         3X1
%           sigma_bar(t)      3X3
function [mu_bar,sigma_bar] = predict(mu,sigma,u,R)
x_r = mu(1:3);
sigma_r = sigma(1:3,1:3);
Gt = [1, 0, -u(2);0,1,u(1);0,0,1];
mu_bar = [(x_r+u),mu(4,end)];
sigma_bar = sigma;
sigma_bar(1:3,1:3) = Gt * sigma_r * (Gt') + R;
end