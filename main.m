clear 
close all
clc
%% run EKF
load map;

runlocalization_track('simout2.txt',M,1,1,1,3);
