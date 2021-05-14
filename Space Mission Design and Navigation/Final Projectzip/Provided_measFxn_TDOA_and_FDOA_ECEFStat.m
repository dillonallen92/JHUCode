%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified:  11-13-2020
% Copyright (c) 2020 Ryan Mitch.  All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code implements an TDOA & FDOA measurement. This follows the course
% reference on TDOA and FDOA.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) x - the [nx x 1] state column vector to use in measurement
%   computation. It corresponds to the time t.  The state is assumed to
%   have the following entries:
%      x = [r_x, r_y, r_z, v_x, v_y, v_z]';
% 
% 2) w - the [nz x 1] measurement noise column vector.  nz = 2.
% 
% 3) t - the scalar time at which the measurement is to be computed.  
% 
% 4) opts - a parameter and general-purpose utility structure.  Two
%   possible inputs are: 
%   1] derFlag - a flag that tells the model to compute the measurement
%       derivatives (Jacobian).  This could be numerical or analytic.  
%       If numeric, then call this same function, but set the flag = 0 to
%       compute each step direction for the state vector.
%       This can be defined outside the function as:  opts.derFlag = ???;
%   2] tFrameAlign - the time that the rotation from ECI to ECEF is the
%       identity matrix.
%       This can be defined outside the function as:  opts.t = ???;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) h - the [nz x 1] column vector of measurements, nz = 2.
% 
% 2) H - the [nz x nx] Jacobian matrix of measurement partial derivatives.
%   The entries are H(j,i) = dh(j)/dx(i).  If you compute numeric
%   derivatives, then finite step methods will work:
%       a) epsilon = 1e-7;  %or similar small number
%       b) compute h (1 function call)
%       c) loop:
%         d) compute x_i = x;
%         e) x(i) = x(i) + epsilon
%         f) compute h_i (call the function again)
%         g) H(:,i) = (h_i – h)/epsilon
% 
% 3) ID - the measurement sensor and number identity.  For now, this can be
%   a dummy value (e.g., 0).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) The measurement noise enters through linear addtion to the
%   measurements.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,H,ID] = Provided_measFxn_TDOA_and_FDOA_ECEFStat(x,w,t,opts)

% Sanity Checks:
nx = size(x,1);
if (nx==1)
    error('x should be a column vector, not a row vector')
end
nz = size(w,1);
if (nz==1) && (size(w,2)~=1)
    error('w should be a column vector, not a row vector')
end

% Default:
ID = 0;

% Measurement equations here:
c     = 299792458;
f0_Hz = opts.f0_Hz;
%
rVec               = x(1:3,1);
rdotVec            = x(4:6,1);
r1_ECEFVec         = opts.r1_ECEFVec;
r2_ECEFVec         = opts.r2_ECEFVec;
dtDay_sidereal_sec = 86164.098903691;
OmegaEarth_radps   = (2*pi)/dtDay_sidereal_sec;
% Station 1:
[r1_IVec]    = Provided_ECI2ECEF_simple(r1_ECEFVec, t, opts.tFrameAlign);
r1_IdotVec   = cross([0;0;OmegaEarth_radps],r1_IVec);
% Station 2:
[r2_IVec]    = Provided_ECI2ECEF_simple(r2_ECEFVec, t, opts.tFrameAlign);
r2_IdotVec   = cross([0;0;OmegaEarth_radps],r2_IVec);
%
% Check elevations:
[elevation1_deg] = computeElevation(rVec,r1_IVec);
[elevation2_deg] = computeElevation(rVec,r2_IVec);
if ~isfield(opts,'ElevationMask_deg')
    ElevationMask_deg = 0;
else
    ElevationMask_deg = opts.ElevationMask_deg;
end
if (elevation1_deg < ElevationMask_deg) || (elevation2_deg < ElevationMask_deg)
    [h,H,ID] = deal([]);    
    return
end

% Compute rho (unit vector):
%1
dr1Vec   = rVec - r1_IVec;
rho1     = sqrt( dr1Vec'*dr1Vec );
rho1Hat  = dr1Vec/rho1;
dv1      = rdotVec - r1_IdotVec;
%2
dr2Vec   = rVec - r2_IVec;
rho2     = sqrt( dr2Vec'*dr2Vec );
rho2Hat  = dr2Vec/rho2;
dv2      = rdotVec - r2_IdotVec;
%
% Eq:
h(1,1) = (1/c)*( norm(dr1Vec) - norm(dr2Vec) );
h(2,1) = (f0_Hz/c)*( (rho1Hat' * dv1) - (rho2Hat' * dv2) );
% Plus noise:
if isempty(w)
    w = 0;
end
h = h + w;

% Jacobian (if needed) here:
if (opts.derFlag == 1)
    % Let's just go numeric here.
    H = [];
    opts.derFlag = 0;
    for ii = 1:6
        xVec_ii = x;
        xVec_ii(ii) = xVec_ii(ii)*(1.0001);
        [hVec_ii,~,~] = Provided_measFxn_TDOA_and_FDOA_ECEFStat(xVec_ii,w,t,opts);
        H(:,ii) = (hVec_ii-h)./(xVec_ii(ii)-x(ii));
    end
else
    % Jacobian not needed.  Save computation time and return.
    H = [];
end