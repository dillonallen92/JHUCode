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
% example problems only- there is not system that produces this kind of
% measurement (GPS is the closes, today).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xVec - The [n x 1] column vector of the state vector.  The first 3 
%   states xVec(1:3) must be the satellite position (frame 
%   agnostic, but must be same as r_IVec).  The other states are ignored.
%  
% 2) w - The [1 x 1] scalar measurement noise for range measurement and any
%   unmodeled propagation delays and errors.  This must be generated 
%   externally to this function.  This is epsilon in book.
% 
% 3) opts.r_IVec - The [3 x 1] column vector of the instrument position 
%   (frame agnostic, but must be same as rVec).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) h - The [3 x 1] position matches.  The units are in the same units as 
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
function [h,HMat,id] = Provided_measFxnRange_ECI_rStatECI(xVec,w,t,opts)
% Parse inputs:
rVec   = xVec(1:3,1);
n      = size(xVec,1);
r_IVec = opts.r_IVec;
% Compute rho:
dr     = rVec - r_IVec;
rho    = sqrt( dr'*dr );
% Compute H:
HMat   = [ (rVec(1)-r_IVec(1))/rho, (rVec(2)-r_IVec(2))/rho,...
    (rVec(3)-r_IVec(3))/rho, zeros(1,(n-3))];

% If truth-model, then add noise:
h = rho;
if ~isempty(w)
    h = h + w;
end

% ID: (not needed yet)
id = [];

return
%% Unit Test:
clear all
close all
clc

rVec = randn(3,1)+1000;
r_IVec = randn(3,1)+(1000*randn(3,1));
xVec = [rVec;randn(3,1)];
n = length(xVec);
w = [];
opts.r_IVec = r_IVec;
[h,H,id] = measFxnRange(xVec,w,opts);
% Check h:
h_true = norm(rVec-r_IVec);
dh = h - h_true;
dh_perc = (dh./h_true)*100
% Check H:
delta = 0.0000001;
for ii = 1:n
    x_ii = xVec;
    x_ii(ii) = x_ii(ii)*(1+delta);
    [h_ii,~,~] = measFxnRange(x_ii,[],opts);
    HNum(1,ii) = (h_ii-h)/(x_ii(ii) - xVec(ii));
end
dH = HNum-H;
dH_perc = (dH./HNum)*100




