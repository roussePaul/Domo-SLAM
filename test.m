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

%% Line extraction

close all
clc

x= linspace(0,1,10);
y=1-x;
X = [x',y'];
S = [0.1,0.01;0.01,0.1]/2;
X = randn(size(X))*S;

plot(X(:,1),X(:,2),'.')


r = linspace(0,5,100);
t = linspace(-pi,pi,100);

[r,t,phi] = lineExtraction(X,S,r,t)

y = (r-cos(t)*x)/sin(t);

hold on
plot(x,y,'r')

%%
clc
clear
close all
N=80;
Nnoise=10



x = linspace(0,0.5,N);
t = -pi;
r=1;
y = (r-cos(t)*x)/sin(t);

X1 = [x',y'];
X1 = X1+0.000005*randn(size(X1));

y = 2-2*x;
X2 = [x',y'];
X2 = X2+0.05*randn(size(X2));

X3 = randn(Nnoise,2)-1;

X = [X1];

tic
[Z,C] = lineExtraction(X,1);
toc
Z(1,:)
C(:,:,1)
x = linspace(-1,1,100);


figure
hold on;
axis equal;

plot(X(:,1),X(:,2),'.');

for i=1:2
    r = Z(i,1);
    t = Z(i,2);

    norm(eig(C(:,:,i)))
    plot(x,(r-x*cos(t))/sin(t));
end

%%

clc
clear
close all
N=80;
Nnoise=10

K=100;
angleIn = linspace(pi/4-0.3,pi/4+0.3,K);
angleOut = angleIn;

for l=1:K

    x = linspace(0,0.5,N);
    t = angleIn(l);
    y = (r-cos(t)*x)/sin(t);

    X1 = [x',y'];
    X1 = X1+0.000005*randn(size(X1));

    y = 2-2*x;
    X2 = [x',y'];
    X2 = X2+0.05*randn(size(X2));

    X3 = randn(Nnoise,2)-1;

    X = [X1];

    tic
    [Z,C] = lineExtraction(X,1);
    toc
    Z(1,:)
    C(:,:,1)
    x = linspace(-1,1,100);

% 
%     figure
%     hold on;
% 
%     plot(X(:,1),X(:,2),'.');
% 
%     for i=1:2
%         r = Z(i,1);
%         t = Z(i,2);
% 
%         norm(eig(C(:,:,i)))
%         plot(x,(r-x*cos(t))/sin(t));
%     end

   angleOut(l) =  Z(1,2);
end

plot(angleIn,angleOut);
p = polyfit(angleIn,angleOut,1)
