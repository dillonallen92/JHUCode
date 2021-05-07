%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified:  5-3-2020
% Copyright (c) 2020 Ryan Mitch.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements a 3D position match measurement.  This is meant for
% example problems only.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  The first 3 
%   states xVec(1:3) must be the satellite position (assumed in inertial).
%   The other states are ignored.
%  
% 2) w - The [1 x 1] scalar measurement noise for range measurement and any
%   unmodeled propagation delays and errors.  This must be generated 
%   externally to this function. 
% 
% 3) opts.r1_ECEFVec - The [3 x 1] column vector of the instrument position 
%   in an Earth-Fixed Frame.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) h - The [1 x 1] range measurement.  The units are in the same units as 
%   r and r_I.
%
% 2) HMat - The [1 x n] matrix of the position measurement partial derivatives 
%   with respect to the states.  
% 
% 3) id - the identification number of the measurement.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Many- this is just for example problems.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,HMat,id] = Provided_measFxnRange_ECI(xVec,w,t,opts)
% Parse inputs:
rVec = xVec(1:3,1);
n    = size(xVec,1);

% Ground station:
r1_ECEFVec  = opts.r1_ECEFVec;
[r_IVec]    = ECI2ECEF_simple(r1_ECEFVec, t, opts.tFrameAlign);

% Compute rho:
dr  = rVec - r_IVec;
rho = sqrt( dr'*dr );
% Compute H: (could speed up by flagging this into a toggleable
% calculation)
HMat = [ (rVec(1)-r_IVec(1))/rho, (rVec(2)-r_IVec(2))/rho,...
    (rVec(3)-r_IVec(3))/rho, zeros(1,(n-3))];

% If truth-model, then add noise:
h = rho;
if ~isempty(w)
    h = h + w;
end

% ID: (not needed yet)
id = [];

return