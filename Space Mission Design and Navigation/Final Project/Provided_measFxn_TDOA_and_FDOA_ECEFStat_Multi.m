%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified:  12-4-2020
% Copyright (c) 2020 Ryan Mitch.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This measurement model implements an ideal TDOA/FDOA measurement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  Inertial frame.
%  
% 2) w - The [(2*nstations) x 1] measurement noise for TDOA/FDOA measurement
%   and any unmodeled propagation delays and errors.  This must be generated 
%   externally to this function. 
% 
% 3) opts - a cell array that has a standard opts input for 
%   Provided_measFxn_TDOA_and_FDOA_ECEFStat.  
%   These should be stored as:  
%       opts{1} = opts for station pair 1
%       opts{2} = opts for station pair 2
%       ...
%   These structures need to have whatever 
%   Provided_measFxn_TDOA_and_FDOA_ECEFStat ask for.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) h - The [(2*nstation) x 1] TDOA & FDOA measurements.  
%
% 2) HMat - The [(2*nstation) x n] matrix of the position measurement partial 
% derivatives with respect to the states.  
% 
% 3) ID - the identification numbers of the measurements.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) if the satellite is below the elevation of either station, then the
% measurement returned is empty.  In other words:  h=[], H=[], ID=[]. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,H,ID] = Provided_measFxn_TDOA_and_FDOA_ECEFStat_Multi(x,w,t,opts)

% Each station gives 2 measurement:
nmeasPerStat = 2;
[h,H,ID] = deal([]);
% Loop over stations:
for ii = 1:numel(opts)
    opts_ii = opts{ii};
    % Grab noise:
    if isempty(w) || (sum(w==0)==numel(w))
        w_ii = 0;
    else
        w_ii = w([1:nmeasPerStat] + (nmeasPerStat*(ii-1)));
    end
    % Evaluate single measurement model:
    [h_ii,H_ii,~] = Provided_measFxn_TDOA_and_FDOA_ECEFStat(x,w_ii,t,opts_ii);
    h = [h; h_ii];
    H = [H; H_ii];
    % Store Meas ID, if valid:
    if isempty(h_ii)
        continue
    else
        ID = [ID; ([1:nmeasPerStat] + ((ii-1)*nmeasPerStat))'];
    end
end
