clear
close all


%% Dataset generation
N=100;

f=1/10;

% simple line
S1 = f*[0.7,0.3;0.3,0.7]/50;
rS1 = sqrtm(S1);
x=linspace(0,1,N);
y=2-3*x;
X = [x',y'];
X = X+randn(size(X))*rS1;
data1 = X;

% static position
S2 = f*[0.7,0.3;0.3,0.7]/5;
rS2 = sqrtm(S2);
x=linspace(0,1,N)*0;
y=2-3*x;
X = [x',y'];
X = X+randn(size(X))*rS2;
data2 = X;

% separate point dataset
S3 = f*[0.7,0.3;0.3,0.7]/20;
rS3 = sqrtm(S3);
p1=[0,0];
p2=[1,1];
X1 = repmat(p1,N,1)+randn(N,2)*rS3;
X2 = repmat(p2,N,1)+randn(N,2)*rS3;
X = [X1;X2];
data3 = X;

% separate line dataset
S4 = f*[0.7,0.3;0.3,0.7]/50;
rS4 = sqrtm(S4);
x1 = linspace(0,1,N);
y1 = 1-x1;
x2 = linspace(-1,0,N);
y2 = 2 + 2*x2;
X = [x1',y1';x2',y2'];
X = X + randn(size(X))*rS4;
data4 = X;


%% Plot
figure
hold on
plot(data1(:,1),data1(:,2),'r.');
plot(data2(:,1),data2(:,2),'g.');
plot(data3(:,1),data3(:,2),'b.');
plot(data4(:,1),data4(:,2),'m.');

%% eigenvalue calculation
d_eig = @(data) norm(eig(cov(data)));
s_eig = @(s) norm(eig(s));
% d_eig(data1)
% s_eig(S1)
% d_eig(data2)
% s_eig(S2)
% d_eig(data3)
% s_eig(S3)

%% Algorithm verification


data = data3;
S = S3;

tic
[Z,C,class] = lineExtraction(data,S,0.3);
toc



figure
hold on
gscatter(data(:,1),data(:,2),class);

N = size(Z,1);
if N==0
    disp('No features detected!');
end

for i=1:N
    x = linspace(min(data(:,1)),max(data(:,1)),2);
    r = Z(i,1);
    t = Z(i,2);
    y = (r-cos(t)*x)/sin(t);
    plot(x,y);
end
