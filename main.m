clear 
close all
clc


%% Options

opt = struct();
opt = setfield(opt,'showEstimate',1);
opt = setfield(opt,'showOdometry',1);
opt = setfield(opt,'showTrue',1);

%% run EKF
load map;

runlocalization_track('simout3.txt',M,opt);