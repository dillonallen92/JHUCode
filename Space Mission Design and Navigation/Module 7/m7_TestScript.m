clear, clc, close all;
load('HW1_Input.mat');

% This script is the test script that will be used to run all the functions
% needed in the homework with the data provided from Blackboard. 

%% Problem 1 (Nav_Sol)

opts.derFlag = derFlag;

[zOrh1, H1, ID] = measFxn_NavSol(x_Fixed, wNavSol, t, opts);


%% Problem 2 (Ideal Range)
x2vec = [r1_IVec; 0;0;0] - x_Inertial;
[zOrh2, H2, ID] = measFxn_IdealRange(x2vec, wIdRng, t, opts);

%% Problem 3 (Ideal Range Rate)
x3vec = [r1_IVec; r1_IdotVec] - x_Inertial;
[zOrh3, H3, ID] = measFxn_IdealRangeRate(x3vec, wIdRngR8, t, opts);

%% Problem 4 (Integrated Carrier Phase)
x4vec = [r1_IVec; r1_IdotVec] - x_Inertial;
x4 = [x4vec; ICP_bias];


%% Problem 5 (Azimuth / Elevation)

%% Problem 6 (TDOA / FDOA)