clear
close all

%% Load data

load data1
% laser : data coming from the laser sensor (361 beams on 180 degrees)
% time : time step oof the measure
% x,y,t : true position of the robot
% v,w : velocity of the robot

load map1
% describe each segment


%% Extract some beams froms the laser sensor
nLaserBeam = 2;

nLaserCSV = size(laser,2);

index = 1:nLaserBeam;
index = floor((index-1)/(nLaserBeam-1)*(nLaserCSV-1))+1;
laser = laser(:,index);

%% Definition of parameters

nbrMeasures = size(laser,1);
nbrSensors = size(laser,2);
nbrPoints = size(laser,1)*size(laser,2);
points = ones(nbrPoints,2);
robotPos = [x,y,t];


%% ideal scenario

j=1;
for i=1:nbrMeasures
    nj=j+nbrSensors-1;
    points(j:nj,:) = getPosFromScan(robotPos(i,:),laser(i,:));
    j=nj;
end

figure
hold on
plot(points(:,1),points(:,2),'.');
plot(x,y);
axis equal;

%% Plot map
N = size(M,3);
for i=1:N
    plot(M(:,1,i),M(:,2,i),'g');
end

%% Compute basic association
N = size(points,1);
class = zeros(N,1);
for i=1:N
    class(i) = associate(M,points(i,:));
end


gscatter(points(:,1),points(:,2),class);
