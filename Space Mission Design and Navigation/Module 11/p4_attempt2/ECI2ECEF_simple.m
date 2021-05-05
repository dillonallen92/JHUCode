%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 12-23-2019
% Copyright (c) 2019 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements a simple ECI to ECEF transformation via a rotation 
% about the z axis.  This is meant to get around the difficulties of
% high-precision frame transformation, since they do not add to the course.
% The frames are assumed to be aligned initially at time tFrameAlign.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) rVec_ECI - [3 x 1] vector of position in the ECI frame
% 
% 2) t - the current time, in seconds.
% 
% 3) tFrameAlign - the time when the ECI/ECEF frames are algined.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) rVec_ECEF - the [3 x 1] vector of position in the ECEF frame.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) frame transformation is a simple function of time and is a rotation
% about the z axis.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) see the NGA appendix on transforming from ECI to ECEF.
% 
% 2) Yes, I checked the direction of rotation.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rVec_ECEF] = ECI2ECEF_simple(rVec_ECI,...
    t, tFrameAlign, varargin)

omega_radps = (2*pi)/86400;
dt_sec = t - tFrameAlign;
% Position:
%       r_ECI = T * r_ECEF
[Rz] = rot(3, (omega_radps*dt_sec));
rVec_ECEF = Rz*rVec_ECI; 

% if ~isempty(varargin)
if 0
    vVec_ECI = varargin{1};
    
    % Velocity:
    %       rdot_ECI = (Tdot*r_ECEF) + (T*rdot_ECEF)
    %
    % % vVec_ECI = (vVec_ECEF) + (omegaVec_radps x rVec_ECEF)
    % % so:  vVec_ECEF = vVec_ECI - (omegaVec_radps x rVec_ECEF)
    %
    %       rdot_ECEF = T' * (rdot_ECI - (Tdot*r_ECEF))
    % omegaVec_radps = [0; 0; omega_radps]; % ??? Not checked.
    % vVec_ECEF = (Rz') * (vVec_ECI - (cross(omegaVec_radps,rVec_ECEF)));
    
    T = dt_sec;
    W = omega_radps;
    dRdt_analytic = [ -W*sin(T*W),  W*cos(T*W), 0;
        -W*cos(T*W), -W*sin(T*W), 0;
        0,           0, 0];
    vVec_ECEF =  (Rz') * (vVec_ECI - (dRdt_analytic*rVec_ECEF))
    warning('Check velocity computation.')
    varargout{1} = vVec_ECEF;
    
    % rVec_ECI = [10e3;0;0]
    % vVec_ECI = [0;8;0]
    %
    % rVec_ECEF = [8192.5; 5734.4; 19]
    % % vECEF = [-4.17; 5.96; 0]
    %
    % omegaVec_radps = [0; 0; (2*pi)/(86400)]; % ??? Not checked.
    % vVec_ECEF = vVec_ECI - (cross(omegaVec_radps,rVec_ECEF))
    % warning('Not checked')
    
    % angle_rad_Vec = linspace(0,2*pi,1e3);
    % angle_rad_Vec = linspace(5.666,5.68,1e3);% best = 5.6725
    % rVec_ECEF_STK = [8192.5; 5734.4; 19];
    % for c = 1:1e3
    %     [Rz] = rot(3, angle_rad_Vec(c));
    %     rVec_ECEF = Rz*rVec_ECI; % Direction?  CHECK IN STK!
    %     J(c) = norm(rVec_ECEF-rVec_ECEF_STK);
    % end
    % so....
    % tFrameAlign = 0;
    % t = angle_rad_Vec(464)/omega_radps = 7.800231593446418e+04
    
    
    
    
    % t1 =  dt_sec;
    % t2 = (dt_sec + 1.0000);
    % [Rz1] = rot(3, t1*(2*pi)/(86400));
    % [Rz2] = rot(3, t2*(2*pi)/(86400));
    % dRdt = (Rz2-Rz1)/(t2-t1)
    % (Rz1')*(vVec_ECI - (dRdt*rVec_ECEF))
    
    
    
    % Ca =cos(T*W)
    % Sa = sin(T*W)
    % Sa =sin(T*W)
    % R = [Ca,Sa,0;
    %         -Sa,Ca,0;
    %         0,0,1]
    %         R =
    % [  cos(T*W), sin(T*W), 0]
    % [ -sin(T*W), cos(T*W), 0]
    % [         0,        0, 1]
end

%% Check rotation over the span of 1 hour vs. STK:
% t0:
% rVec_ECI = [10e3;0;0]
% vVec_ECI = [0;8;0]
% rVec_ECEF = [8192.5; 5734.4; 19]
% vECEF = [-4.17; 5.96; 0]
% t0 = t1 + 1 hour
% rVec_ECI = [10e3;0;0]
% vVec_ECI = [0;8;0]
% rVec_ECEF = [9.4e3; 3.4e3; 0.019]
% vVec_ECEF = [-2.5; 6.8; 0]
%
%
% [rVec_ECEF] = ECI2ECEF_simple(rVec_ECI, t, tFrameAlign)
% rVec_ECEF =   1.0e+03 *[8.1925;   5.7344;   0];
% [rVec_ECEF] = ECI2ECEF_simple(rVec_ECI, t+3600, tFrameAlign)
% rVec_ECEF =   1.0e+03 *[9.3975;   3.4186;   0];


