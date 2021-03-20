%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dillon Allen
% Date Modified: 3/19/2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code implements an <???> measurement. This implementation follows 
% equations out of Tapley, Schutz, and Born Chapter <???>(2004 edition) and 
% the course notes for Module 7.
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
% 2) w - the [nz x 1] measurement noise column vector.
% 
% 3) t - the scalar time at which the measurement is to be computed.  This
%   may or may not be used.  In measurements that require frame rotations 
%   (e.g., ECEF to ECI) this time is critical.
% 
% 4) opts - a parameter and general-purpose utility structure.  Two
%   possible inputs are: 
%   1] derFlag - a flag that tells the model to compute the measurement
%       derivatives (Jacobian).  This could be numerical or analytic.  
%       If numeric, then call this same function, but set the flag = 0 to
%       compute each step direction for the state vector.
%       This can be defined outside the function as:  opts.derFlag = ???;
%   2] t0 - the time that the rotation from ECI to ECEF is the identity
%       matrix.
%       This can be defined outside the function as:  opts.t = ???;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) zOrh - the [nz x 1] column vector of measurements (simulated, or
%   predicted).  z uses a non-zero w and h uses a zero w (this is something
%   you determine outside of this function).
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
function [zOrh,H,ID] = measFxn_IdealRangeRate(x,w,t,opts)

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
%  zOrh = ???


% Jacobian (if needed) here:
if (opts.derFlag == 1)
    % Compute Jacobian:
    % H = ???
else
    % Jacobian not needed.  Save computation time and return.
    H = [];
end