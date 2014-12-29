clear 
close all
clc

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


%% Options

opt = struct();
opt = setfield(opt,'showEstimate',1);
opt = setfield(opt,'showEstimateCov',1);
opt = setfield(opt,'showOdometry',0);
opt = setfield(opt,'showTrue',1);
opt = setfield(opt,'showStep',500);
opt = setfield(opt,'verbose',3);
opt = setfield(opt,'maxStep',inf);
opt = setfield(opt,'trueMap',M);
opt = setfield(opt,'showTrueMap',1);


%% run EKF
load map;

runlocalization_track('simout_5_sensors_noisy.txt',M,opt);

