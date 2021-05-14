%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified:  12-4-2020
% Copyright (c) 2020 Ryan Mitch.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This measurement model implements an ideal range measurement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  Inertial frame.
%  
% 2) w - The [nstations x 1] measurement noise for range measurement and any
%   unmodeled propagation delays and errors.  This must be generated 
%   externally to this function. 
% 
% 3) opts - a cell array that has a standard opts input for 
%   Provided_measFxnRange_ECI.  
%   These should be stored as:  
%       opts{1} = opts for station 1
%       opts{2} = opts for station 2
%       ...
%   These structures need to have whatever Provided_measFxnRange_ECI asks
%   for, namely r1_ECEFVec - the [3 x 1] column vector of the instrument 
%   position in an Earth-Fixed Frame, and the time of frame align.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) h - The [nstation x 1] range measurement.  
%
% 2) HMat - The [nstation x n] matrix of the position measurement partial 
% derivatives with respect to the states.  
% 
% 3) ID - the identification numbers of the measurements.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) if the satellite is below the elevation, then the measurement returned
% is empty.  In other words:  h=[], H=[], ID=[]. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,H,ID] = Provided_measFxnRange_ECI_multi(xVec,w,t,opts)

% Each station gives 1 measurement:
nmeasPerStat = 1;
[h,H,ID] = deal([]);
% Loop over stations:
for ii = 1:numel(opts)
    opts_ii = opts{ii};
    % Grab noise:
    if isempty(w) || (sum(w==0)==numel(w))
        w_ii = 0;
    else
        w_ii = w(ii);
    end
    % Evaluate single measurement model:
    [h_ii,H_ii,~] = Provided_measFxnRange_ECI(xVec,w_ii,t,opts_ii);
    h = [h; h_ii];
    H = [H; H_ii];
    % Store Meas ID, if valid:
    if isempty(h_ii)
        continue
    else
        ID = [ID; ([1:nmeasPerStat] + ((ii-1)*nmeasPerStat))'];
    end
end