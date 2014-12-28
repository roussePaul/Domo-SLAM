%% draw state
clc
close all
clear

mu = [0;0;pi*2/3;3.4;pi/3];
sensorPose = [0.1;0];

axis equal

drawState(mu,sensorPose,[-1,1,-1,1]);

drawMeasure(mu,sensorPose, [0.4,pi;0.5,pi/2;0.3,pi/4]);

%% observation_model
clc
close all
clear

mu = [1;0;pi/3;0.3;0.3];
sensorPose = [0.2;0];
a = [-1,2,-1,2];
drawState(mu,sensorPose,a);

t = -pi/2;
z = observation_model(mu,1,t)

drawMeasure(mu,sensorPose,[z,t])

%% jacobian observation model
clc

dx = 1e-8;

mu = randn(27,1);

z = observation_model(mu,1,pi/2);

H = jacobian_observation_model(mu,1,z,1,pi/2);
i=0;
H2 = H;
for i=0:4
    mu_d = (mu+dx*[zeros(i,1);1;zeros(26-i,1)]);
    zdx = observation_model(mu_d,1,pi/2);
    H2(i+1) = -(z-zdx)/dx;
end

H(1:5)
H2(1:5)
err = (H-H2)./H;
err(1:5)
