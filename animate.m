% Animation of robot
clear
close all
clc

%% load data
load animate;

%% animate
hold on;
f = drawFeature(mu,[0.2;0],[-4 48 -33 9]);
r=[];
m=[];
for i=1:N
    delete(r);
    delete(m);
    mu = [x(i);y(i);t(i);mu(4:end)];
    r = drawRobot(mu,[0.2;0],[-4 48 -33 9]);
    m = drawMeasure(mu ,[0.2;0], [laser(i,1),-pi/2;laser(i,2),pi/2]);
    h = [f r m];
    drawnow;
end
