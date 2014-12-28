clear 
close all
clc


%% Options

opt = struct();
opt = setfield(opt,'showEstimate',0);
opt = setfield(opt,'showEstimateCov',1);
opt = setfield(opt,'showOdometry',0);
opt = setfield(opt,'showTrue',0);
opt = setfield(opt,'verbose',1);
opt = setfield(opt,'maxStep',5000);

%% run EKF
load map;

runlocalization_track('simout.txt',M,opt);