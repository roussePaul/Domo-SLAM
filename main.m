clear 
close all
clc
%% run EKF
load map;

runlocalization_track('simout3.txt',M,1,0,0,3);
