% Animation of robot
clear
close all
clc

%% load data
load animate;

%% animate
hold on;
hf = drawFeature(mu,[0.2;0],[-4 48 -33 9]);
hr=[];
hm=[];
for i=1:N
    delete(hr);
    delete(hm);
    mu = [x(i);y(i);t(i);mu(4:end)];
    hr = drawRobot(mu,[0.2;0],[-4 48 -33 9]);
    hm = drawMeasure(mu ,[0.2;0], [laser(i,1),-pi/2;laser(i,2),pi/2]);
    h = [hf hr hm];
    drawnow;
end
