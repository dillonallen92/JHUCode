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
function [zOrh,H,ID] = measFxn_TDOA_FDOA(x,w,t,opts)

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

% Constant
c = 1e8;

r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];

r1 = [opts.r1(1); opts.r1(2); opts.r1(3)];
v1 = [opts.r1(4); opts.r1(5); opts.r1(6)];

r2 = [opts.r2(1); opts.r2(2); opts.r2(3)];
v2 = [opts.r2(4); opts.r2(5); opts.r2(6)];

freq = opts.f0;

rho1 = (1/norm(r-r1)) * (r-r1);
rho2 = (1/norm(r-r2)) * (r-r2);
% Measurement equations here:
%  zOrh = ???
zOrh(1) = (1/c)*(norm(r-r1) - norm(r-r2)) + w(1);
zOrh(2) = (freq/c)*(dot(rho1,(v-v1)) - dot(rho2, (v-v2))) + w(2);
 
% Jacobian (if needed) here:
if (opts.derFlag == 1)
    % Compute Jacobian:
    % H = ???
    syms rx ry rz vx vy vz r1x r1y r1z v1x v1y v1z r2x r2y r2z v2x v2y v2z;
    r = [rx; ry; rz];
    v = [vx; vy; vz];
    r1 = [r1x; r1y; r1z];
    v1 = [v1x; v1y; v1z];
    r2 = [r2x; r2y; r2z];
    v2 = [v2x; v2y; v2z];
    rho1 = (1/norm(r-r1)) * (r-r1);
    rho2 = (1/norm(r-r2)) * (r-r2);
    zMat = [(1/c)*(norm(r-r1) - norm(r-r2)); (freq/c)*(dot(rho1,(v-v1)) - dot(rho2, (v-v2)))];
    variableVec = [rx ry rz vx vy vz r1x r1y r1z v1x v1y v1z r2x r2y r2z v2x v2y v2z];
    numericalVec = [x(1) x(2) x(3) x(4) x(5) x(6) opts.r1(1) opts.r1(2) opts.r1(3) opts.r1(4) opts.r1(5) opts.r1(6) ...
                    opts.r2(1) opts.r2(2) opts.r2(3) opts.r2(4) opts.r2(5) opts.r2(6)];
    H = eval(subs(jacobian(zMat, [rx ry rz vx vy vz]), variableVec, numericalVec));
    
else
    % Jacobian not needed.  Save computation time and return.
    H = [];
end