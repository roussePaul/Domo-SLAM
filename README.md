Domo-SLAM
=========

SLAM implementation with indoor environement geometric constraints.


Basic use
=======
Just execute the main.m file in the project directory.


Create Dataset
===========
Use DatasetGeneration.m

One dataset have already been created, they are in the simout/ directory.


Configure the SLAM algorithm
=====================

main.c
Here you can configure the output (trajectories plot, error, covariance plots...)

init.m
You can change the configuration of the initial state and covariance matrix.

addFeatures.m
Configuration of the line extraction algorithm.