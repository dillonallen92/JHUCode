clear, clc, close all;
load('HW1_Input.mat');

% This script is the test script that will be used to run all the functions
% needed in the homework with the data provided from Blackboard. 

%% Problem 1 (Nav_Sol)

opts.derFlag = derFlag;

[zOrh, H, ID] = measFxn_NavSol(x_Fixed, wNavSol, t, opts);


%% Problem 2 (Ideal Range)

%% Problem 3 (Ideal Range Rate)

%% Problem 4 (Integrated Carrier Phase)

%% Problem 5 (Azimuth / Elevation)

%% Problem 6 (TDOA / FDOA)