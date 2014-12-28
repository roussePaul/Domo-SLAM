% function [mu,sigma,R,Q,Lambda_M] = init()
% This function initializes the parameters of the filter.
% Outputs:
%			mu(0):			3X1
%			sigma(0):		3X3
%			R:				3X3
%			Q:				2X2
function [mu,sigma,R,Q,Lambda_M] = init()
load mu_init mu;
sigma = 1e-4*diag(ones(size(mu))); % initial covariance matrix
delta_m = 0.999;
Lambda_M = chi2inv(delta_m,2);

R = diag([0.1^2, 0.1^2, 0.01^2])/100;
Q = diag([0.1^2])*1000;

end