%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 10-20-2019
% Copyright (c) 2019 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements a conversion from the Earth Centered Earth Fixed
% (ECEF) frame to the Tangent East North Vertica (ENV) frame.
% 
% Core Algorithm:  
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% A unit test has been included at the bottom.
% 
%  difficulty- easy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rENV_m, TECEF2ENV] = frame_ECEF2TangentENV(rECEF_m, rOriginECEF_m)

% Offset of frame center:
rPrime_m = rECEF_m - rOriginECEF_m;
%
% Rotation matrix for the new origin point:
rho_m = sqrt((rOriginECEF_m(1).^2) + (rOriginECEF_m(2).^2));
lat_rad = atan2(rOriginECEF_m(3), rho_m);
lon_rad = atan2(rOriginECEF_m(2), rOriginECEF_m(1));
%
[Rz] = rot(3, lon_rad);
[Ry] = rot(2, -lat_rad);
%
T_ECEF2ENV_B4Permute = Ry*Rz;
P = [0,1,0; 0,0,1; 1,0,0]; % permutation matrix to reorder to ENV
TECEF2ENV = P*T_ECEF2ENV_B4Permute;
%
%
% Change frame:
rENV_m = TECEF2ENV*(rPrime_m);




return

%% Unit test:
clear all 
clc

% 0 lat,lon
Re_m = 6371e3;
lat_deg = 0;
lon_deg = 0;
r = Re_m * [(cosd(lat_deg)*cosd(lon_deg));  cosd(lat_deg)*sind(lon_deg);  sind(lat_deg)];
% Unchanged?
[rENV_m] = frame_ECEF2TangentENV(r,r)


% 90 lat, 0 lon
Re_m = 6371e3;
lat_deg = 89;
lon_deg = 0;
r = Re_m * [(cosd(lat_deg)*cosd(lon_deg));  cosd(lat_deg)*sind(lon_deg);  sind(lat_deg)];
% 
[rENV_m] = frame_ECEF2TangentENV(r+[1;0;0],r) %-1 in N
[rENV_m] = frame_ECEF2TangentENV(r+[0;1;0],r) %+1 in E
[rENV_m] = frame_ECEF2TangentENV(r+[0;0;1],r) %+1 in V

% 0 lat, 90 lon
Re_m = 6371e3;
lat_deg = 0;
lon_deg = 90;
r = Re_m * [(cosd(lat_deg)*cosd(lon_deg));  cosd(lat_deg)*sind(lon_deg);  sind(lat_deg)];
% 
[rENV_m] = frame_ECEF2TangentENV(r+[1;0;0],r) %-1 in E
[rENV_m] = frame_ECEF2TangentENV(r+[0;1;0],r) %+1 in V
[rENV_m] = frame_ECEF2TangentENV(r+[0;0;1],r) %+1 in N


% 89 lat, 270 lon
Re_m = 6371e3;
lat_deg = 89;
lon_deg = 270;
r = Re_m * [(cosd(lat_deg)*cosd(lon_deg));  cosd(lat_deg)*sind(lon_deg);  sind(lat_deg)];
% 
[rENV_m] = frame_ECEF2TangentENV(r+[1;0;0],r) %+1 in E
[rENV_m] = frame_ECEF2TangentENV(r+[0;1;0],r) %+1 in N
[rENV_m] = frame_ECEF2TangentENV(r+[0;0;1],r) %+1 in V




