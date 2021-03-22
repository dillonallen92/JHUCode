clear, clc, close all;
load('HW1_Input.mat');

% This script is the test script that will be used to run all the functions
% needed in the homework with the data provided from Blackboard. 

%% Problem 1 (Nav_Sol)
opts.derFlag = derFlag;
[zOrh1, H1, ID] = measFxn_NavSol(x_Fixed, wNavSol, t, opts);

%% Problem 2 (Ideal Range)
x2vec =  x_Inertial - [r1_IVec; 0;0;0];
[zOrh2, H2, ID] = measFxn_IdealRange(x2vec, wIdRng, t, opts);

%% Problem 3 (Ideal Range Rate)
x3vec = x_Inertial - [r1_IVec; r1_IdotVec];
[zOrh3, H3, ID] = measFxn_IdealRangeRate(x3vec, wIdRngR8, t, opts);

%% Problem 4 (Integrated Carrier Phase)
x4vec = x_Inertial - [r1_IVec; r1_IdotVec];
x4 = [x4vec; ICP_bias];
[zOrh4, H4, ID] = measFxn_IntegratedCarrierPhase(x4, wIntCarPhs, t, opts);

%% Problem 5 (Azimuth / Elevation)
x5vec = x_Inertial - [r1_IVec; 0; 0; 0];
opts.t = tFrameAlign;
[zOrh5, H5, ID] = measFxn_Azimuth_Elevation(x5vec, wAzEl, t, opts);

%% Problem 6 (TDOA / FDOA)
% Not sure how to input both vectors, so I guess use the opts since its
% the general parameter structure?
x6vec = x_Inertial;
opts.r1 = [r1_IVec; r1_IdotVec];
opts.r2 = [r2_IVec; r2_IdotVec];
opts.f0 = f0_Hz;
opts.t = tFrameAlign;

[zOrh6, H6, ID] = measFxn_TDOA_FDOA(x6vec, wTDOAFDOA, t, opts);

