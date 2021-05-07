function [tout,xOut,Phi,Gu,Gv] = Provided_dynFxnHand_singleClockModel_ZT(x,u,...
    v,tV,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a simple implementation of the Zucca and Tavella clock model, in
% a dynamics model form suitable to filtering or linear covariance
% analysis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   x is the state column vector, 
%   u is the disturbance input (0 unless you want to simulate disturbances)
%   v is the  control input (0 unless you want to simulate contrl voltage changes) 
%   tVec is the time vector [tfrom, tto], and
%   opts is an options input.  
% 
% 
% Single clock dynamics following the paper by Zucca and Tavella:
%   assumes mu_3 is either in the parameter space or is set to 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics:
tau = tV(end) - tV(1);
Phi = [1,tau,(tau^2)/2;...
       0,1,tau;...
       0,0,1];
   
% Disturbance:
Gu = eye(3);  

% Control:
Gv = eye(3);

% Combine Dynamics:
xOut = (Phi * x) + (Gu * u) + (Gv * v);

% Output:
tout = tV(end);


