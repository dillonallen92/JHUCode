%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 03-29-2021
% Copyright (c) 2020 Ryan Mitch/JHUAPL
% Template Modified by: Dillon Allen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes a dynamics model for use in later estimation
% routines.  Requires as output the 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) x - the (nx x 1) state vector at time k
% 2) u - the (nu x 1) disturbance input vector at time k.  It is assumed 
%       to be held constant throughout the integration interval.
% 3) v - the (nv x 1) control input vector at time k.  It is 
%       assumed to be held constant throughout the integration interval.
% 4) tspan- the time span over which to evaluate the dynamics.  This
%       can be a two element vector [t_start, t_end] or a column vector.
%       The units are assumed to be in seconds.
% 5) param - this is a parameter structure that contains all of the other
%       inputs that might be needed by various dynamics models. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1) t - the final time that the state is defined at, in seconds.
%  2) x - the (nx x 1) state vector defined at time tspan(end). 
%  3) Phi - the (nx x nx) state transition matrix.  If it were a linear 
%           system, then the STM would take a state from time t_k and 
%           propagates it to time t_kp1.  It is just an approximation.
%  4) Gamma_u - the (nx x nu) process noise influence matrix that takes a
%           disturbance input vector and defines its influence from 
%           t_k->t_kp1.
%  5) Gamma_v - the (nx x nv) control influence matrix that takes a
%           control input vector and defines its influence from t_k->t_kp1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,Phi,Gamma_u,Gamma_v] = dynamicsFunc(x,u,v,tspan,param)

% Sanity Checks:
nx = size(x,1);
if (nx==1)
    error('x should be a column vector, not a row vector')
end

% Time output:
t = tspan(end);

% Update x:  (IF linear, x = Phi*x + Gamma_u*u + Gamma_v*v)


% Jacobian (if needed) here:
if (param.derFlag == 1)
    % Compute Phi:
    % Phi = ???;
    Phi = @(t,s)[1 t-s (t-s)^2 / 2; 0 1 t-s; 0 0 1];
    % Compute Gamma_u and Gamma_v:
    Gamma_u = [1 0 0; 0 1 0; 0 0 1];    
    Gamma_v = [ 1 0 0; 0 1 0; 0 0 1];    
else
    % Jacobian not needed.  Save computation time and return.
    [Phi,Gamma_u,Gamma_v] = deal([]);
end

% Update x:  (IF linear, x = Phi*x + Gamma_u*u + Gamma_v*v)
x(:,1) = x;
for i = 1 : length(tspan)-1
   x(:,i+1) = Phi(tspan(i+1), tspan(i))*x(:,i) + Gamma_u * u(:,i) + Gamma_v * v;
end

end


