%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 1-10-2021
% Copyright (c) 2021 Ryan Mitch.  All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements an azimuth elevation measurement. This 
% implementation follows equations of Tapley, Schutz, and Born 
% Chapter 2.4.3 (2004 edition).  This assumes the instrument/reference is
% at the surface of the Earth and the satellite state vector is defined in
% the ECI coordinate system.
% 
% Core Algorithm:  
%  Azimuth & Elevation (use atan2 for Az)
%       sin(EL) = z/r
%       sin(Az) = x/r_xy
%       cos(Az) = y/r_xy
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  The first 3 
%   states xVec(1:3) must be the satellite position in the ECI
%   frame.  The other states are ignored.
%  
% 2) wVec - The [nz x 1] measurement noise for an azimuth and elevation
%   measurement and any unmodeled propagation delays and errors.  
%   This must be generated externally to this function.  
% 
% 3) t - the effective time of the measurement.
% 
% 4) opts - a parameter and general-purpose utility structure.  Two
%   possible inputs are: 
%   1] derFlag - a flag that tells the model to compute the measurement
%       derivatives (Jacobian).  This could be numerical or analytic.  
%       If numeric, then call this same function, but set the flag = 0 to
%       compute each step direction for the state vector.
%       This can be defined outside the function as:  opts.derFlag = ???;
%   2] tFrameAlign - the time that the rotation from ECI to ECEF is 
%       the identity matrix.
%   3] r_IVec - the observing station location (ECI frame).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) h - The [2 x 1] vector of the predicted azimuth and then elevation 
%   of the satellite wrt the instrument, given the inputs.  The units 
%   are radians.
%
% 2) HMat - The [1 x n] matrix of the az/el  measurement partial 
%   derivatives with respect to the states.  
% 
% 3) id - the identification number of the measurement.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) instantaneous az/el
% 2) no propagation delays
% 3) no unmodeled errors (note the use of wVec above)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation here follows the TSB book and it is meant to be
% understandable to the reader (for class).  There are more efficient and
% stable implementations of this algorithm.  
% 
% Az and El measurements have singularities/jumps that can cause issues.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,HMat,id] = Provided_measFxnAzEl_ECI_HW11(xVec,wVec,t,opts)

% Ground station:
r_ECIVec         = xVec(1:3);
[r_ECEFVec]      = ECI2ECEF_simple(r_ECIVec, t, opts.tFrameAlign);
rStationECEF_m   = opts.rStationECEF_m;

% Opts must contain the location of the observing station:
[rENV_m, TECEF2ENV] = frame_ECEF2TangentENV(r_ECEFVec,rStationECEF_m);

% Parse inputs:
r_tVec = rENV_m;
n      = size(xVec,1);

% Compute azimuth and elevation
r_t    = norm(r_tVec);
% r_t_xy = norm(r_tVec(1:2));
%
arg_el = r_tVec(3)/r_t;
%
% az_rad = asin(arg_az);% not complete
az_rad = atan2(r_tVec(1),r_tVec(2));
el_rad = asin(arg_el);
% If truth-model, then add noise:
h = [az_rad;...
    el_rad];

% Check edge/acceptable conditions:
if ( isnan(h(1)) || (isnan(h(2))) )
    error('Edge case detected.')
end
% Apply elevation mask if below 0 degrees elevation (cannot be seen!):
if (el_rad<0)
    [h,HMat,id] = deal([]);
    return
else   
    % ID: (not needed yet)
    id = [];    
end
%
% H:
HMat = [];
if isfield(opts,'derFlag')
    if (opts.derFlag==1)
        %
        opts.derFlag = 0;
        % Compute partial derivatives:
        HMat = zeros(2,n);
        % Numerically solve for the derivative (the math is messy):
        delta = 1e-7;
        for ii = 1:3        
            x_ii = xVec;
            x_ii(ii) = x_ii(ii)*(1+delta);
            [h_ii,~,~] = Provided_measFxnAzEl_ECI_HW11(x_ii,[],t,opts);
            try
            HMat(:,ii) = (h_ii-h)/(x_ii(ii) - xVec(ii));
            catch
                keyboard
            end
        end
    end
end

if ~isempty(wVec)
    if (size(wVec,1)==2) && (size(wVec,2)==1)
        h = h + wVec;
    else
        error('Sizing issue on measurement noise.')
    end
end

