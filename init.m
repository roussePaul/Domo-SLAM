% function [mu,sigma,R,Q,Lambda_M] = init()
% This function initializes the parameters of the filter.
% Outputs:
%			mu(0):			3X1
%			sigma(0):		3X3
%			R:				3X3
%			Q:				2X2
function [mu,sigma,R,Q,Lambda_M] = init()
load mu_init;   % load mu_alpha the initial vector for the SLAM algorithm. 
                % If you want to have knowledge about the map, don't change mu_alpha

                
mu = mu_alpha(1:4); % We remove any feature information and we keep the alpha angle of the map

sigma = 1e-4*diag(ones(size(mu))); % initial covariance matrix of the EKF

sigma(1:3,1:3) = 0;     % It is a SLAM algorithm, therefore we assume that the robot is at a already known position at the beginning
sigma(4,4) = 1e-7;      % We set \Sigma_{\alpha} with a good belief, otherwise the EKF diverge!

delta_m = 0.999;
Lambda_M = chi2inv(delta_m,2);


R = diag([0.1^2, 0.1^2, 0.01^2])/100;   % Noise covariance matrix of the process
Q = diag([1^2]);            % Noise covariance matrix of the measure

end