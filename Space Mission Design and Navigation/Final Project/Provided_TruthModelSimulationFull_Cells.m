%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 12-4-2020
% Copyright (c) 2020 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a truth model simulation that generates measurement and process
% noise to create realistic statistics.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) dynFxnHand - the dynamics model function handle.  This is meant to
% generalize the function.  A standard MATLAB function can be called by
% definining a variable for that function:  
%   functionHandleVariable = @functionName
% In this case the function is supposed to follow this format:
%   [tout,xOut,F,Gu,Gv] = dynFxnHandTest3(x,u,v,tV,opts) 
%   where x is the state vector, u is the disturbance input, v is a
%   control input (0 here), tVec is the time vector [tfrom, tto], and
%   opt is an options function.  
%   tout is the output time (should be tV(end)), xOut is the state after
%   propagation, F is the integrated state transistion matrix, Gu is the
%   process noise influence matrix, Gv is the control influence
%   matrix (0 here).
% 
% 2) dynOpts - the dynamics model options.
% 
% 3) measFxnHand - the measurement model function handle.  Similar to
% dynFxnHand, but for the measurements.  In this case it should follow this
% format:
%  [h,H,id] = measFxnHand(x,w,t,opts)
%   where x is the state after propagation, w is the noise vector (should
%   be 0 in estimation, but non-zero in truth-model simulation.
% 
% 4) measOpts - the measurement function handle.
% 
% 5) xTrue0 - the [nx x 1] true state vector at the initial time.
% 
% 6) tVec - The [(k+1) x 1] vector of times to simulate the state for and
%    measurements for.  The first time is the time of the initial state 
%   (tm1) and does not correspond to a measurement time.  All other times 
%   corrspond to a measurement time.
% 
% 7) RMat - The [m x m] matrix of measurement noise covariance matrix. 
% 
% 8) QMats - the [n x n] process noise covariance matrices for each of
% the k times steps.  The process noise (u) is assumed to be 0 mean (E[u]=0) 
% and have Q covariance (E[(u-ubar)(u-ubar')] = Q).  This needs to be 
% selected at the same time as the dynamcis model, as the dynamics model's
% Gu will propagate that process noise.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) zCellHist - The cell array of measurements vs time.
%
% 2) zNoNoiseCellHist - The cell array of measurements vs time, but without
%   corruption by measurement noise.
% 
% 3) xTrueCellHist - The true state time history stored in cells.
%  
% 4) IDsCells - The IDs of the measurements.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) x0 is at t0 (tVec(1))
% 2) RMat is time-invariant
% 3) QMat is time-invariant
% 4) not control (vInfo=[])
% Format of dyn/meas defined in notationConvention.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zCellHist,zNoNoiseCellHist,xTrueCellHist,IDsCells] =...
    Provided_TruthModelSimulationFull_Cells(dynFxnHand,dynOpts,measFxnHand,...
    measOpts,xTrue0,tVec,RMat,QMat)
%
nk = length(tVec)-1;
xTrue = xTrue0;
% Initialize Storage:
IDsCells = cell(nk,1);
xTrueCellHist = cell(nk,1);
xTrueCellHist{1} = xTrue0;
for k = 1:nk
    % Dynamics:
    tk = tVec(k);
    tkp1 = tVec(k+1);
    u = chol(QMat)*randn(size(QMat,1),1);
    [~,xTruePropVec,~,~,~] = feval(dynFxnHand,xTrue,u,[],[tk,tkp1],dynOpts);
    xTrue = xTruePropVec(:,end);
    %
    % Observations:
    wVec = chol(RMat)*randn(size(RMat,1),1);
    % Common failure point:
    try
        [zCellHist{k},~,IDs] = feval(measFxnHand,xTrue,wVec,tkp1,measOpts);
    catch
        keyboard
    end
    [zNoNoiseCellHist{k},~,~] = feval(measFxnHand,xTrue,[],tkp1,measOpts);
    IDsCells{k} = IDs;
    %
    % Storage:
    xTrueCellHist{k+1} = xTrue;
end