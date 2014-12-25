% function H = jacobian_observation_model(mu_bar,M,j,z,i)
% This function is the implementation of the H function
% Inputs:
%           mu_bar(t)   (3+2*N)X1
%           j           1X1 feature index (between 1 and N=nbre of
%           features)
%           z           2Xn
%           psi         1X1 angle of the sensor
% Outputs:  
%           H           1x(3+2*N)
function H = jacobian_observation_model(mu_bar,j,z,psi)
j_rho = 3+2*j-1;
j_theta = j_rho + 1;

N = (size(x,1)-3)/2;    % Nbre of features

rho = mu_bar(j_rho);
theta = mu_bar(j_theta);

phi = x(3);
xr = x(1);
yr = x(2);

r = 1/sin(theta - psi - phi);

H = [
    r*cos(theta),   %dz/dx
    r*sin(theta),   %dz/dy
    -d*sin(phi-theta)*r + cos(theta-psi-phi) * z * r,   %dz/dphi
    zeros(1,2*j-2),
    -r,             %dz/drho_i
    r*(-xr*sin(theta)+yr*cos(theta) + d*sin(phi-theta)) - cos(theta-psi-phi) * z *r,    %dz/dtheta_i
    zeros(1,2*(N-j))
    ];
end
