%% Dataset generation

clear
close all

%% Load data

load data1
% laser : data coming from the laser sensor (361 beams on 180 degrees)
% time : time step oof the measure
% x,y,t : true position of the robot
% v,w : velocity of the robot

%% Map
M = zeros(2,2,12);
i=1;

M(:,:,i) = [-2.76,-6.40;-3.75,2.84]; i=i+1;
M(:,:,i) = [M(2,:,i-1);44.01,8.303]; i=i+1;
M(:,:,i) = [M(2,:,i-1);47.97,-26.82]; i=i+1;
M(:,:,i) = [M(2,:,i-1);25.86,-29.3]; i=i+1;
M(:,:,i) = [M(2,:,i-1);24.76,-19.55]; i=i+1;
M(:,:,i) = [M(2,:,i-1);12.71,-20.89]; i=i+1;
M(:,:,i) = [M(2,:,i-1);13.78,-30.62]; i=i+1;
M(:,:,i) = [M(2,:,i-1);0.1436,-32.18]; i=i+1;
M(:,:,i) = [M(2,:,i-1);-2.647,-7.122]; i=i+1;
M(:,:,i) = [M(2,:,i-1);19.14,-4.681]; i=i+1;
M(:,:,i) = [M(2,:,i-1);19.09,-4.036]; i=i+1;
M(:,:,i) = [M(2,:,i-1);M(1,:,1)];


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
N = size(laser,1);
j=1;
for i=1:N
    nj=j+nbrSensors;
    points(j:(nj-1),:) = getPosFromScan(robotPos(i,:),laser(i,:));
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
    class(i) = associate2(M,points(i,:));
end
gscatter(points(:,1),points(:,2),class,'bgrcmyk');





%% Create file

file = fopen('simout.txt','w');
N = size(laser,1);
for i=1:N
    fprintf(file,'%f %f %f %f %f %f %d %f %f %d %f %f %d\n',...
        time(i), x(i), y(i), t(i), v(i), w(i), 2,...
        pi/2, laser(i,1), class(i), -pi/2,  laser(i,2), class(N+i) );
end

fclose(file);

%% Compute initial state (map and position)
N = size(M,3);
mu = zeros(2*N,1);


for i=1:N
    theta = atan2(M(2,2,i)-M(1,2,i),M(2,1,i)-M(1,1,i))-pi/2;
    theta = mod(theta+pi,2*pi)-pi;
    rho = M(1,1,i)*cos(theta) + M(1,2,i)*sin(theta);
    mu(2*i-1) = rho;
    mu(2*i) = theta;
end

mu  =[x(1);y(1);t(1);mu];

save mu_init mu