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