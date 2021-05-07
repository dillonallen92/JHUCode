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
% example problems only.  The match is performed in the ECEF (state is
% assumed ECI).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  
%  
% 2) w - The [3 x 1]  measurement noise for pos fix measurement and any
%   unmodeled propagation delays and errors.  This must be generated 
%   externally to this function.  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) h - The [3 x 1] position matches.  The units are in the same units as 
%   x.
%
% 2) HMat - The [1 x n] matrix of the position measurement partial 
%   derivatives with respect to the states.  
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
function [horZ,H,id] = Provided_measFxnNavSol_ECI(xVec,w,t,opts)
% Parse inputs:
rVec = xVec(1:3,1);

% Ground station:
[rECEFVec] = ECI2ECEF_simple(rVec, t, opts.tFrameAlign);

% Check altitude:
Re_m        = 6.371e6;
alt_m       = norm(rECEFVec) - Re_m;
altCutOff_m = 20e6;
if (alt_m > altCutOff_m)
    [horZ,H,id] = deal([]);
    return
end

% Measurement:
horZ = rECEFVec;
if ~isempty(w)
    horZ = horZ + w;
end

if (opts.derFlag == 1)
    % Compute Jacobian: 
    H = [];
    opts.derFlag = 0;
    epsilon = 1e-6;
    for ii = 1:6
        xVec_ii = xVec;
        xVec_ii(ii) = xVec_ii(ii)*(1+epsilon);
        [hVec_ii,~,~] = Provided_measFxnNavSol_ECI(xVec_ii,w,t,opts);
        try
        H(:,ii) = (hVec_ii-horZ)./(xVec_ii(ii)-xVec(ii));
        catch
            keyboard
        end
    end
else
    % Jacobian not needed.  Save computation time and return.
    H = [];
end

% ID: (not needed yet)
id = [1:3]';


