%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 10-26-2020
% Copyright (c) 2020 Ryan Mitch/JHUAPL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function executes the dynamics defined in a continuous time domain
% and converts them into the discrete time approximations through a
% difference equation.  The general form of which is:
%
% Continuous time:
% dxdt = (A*x) + (B_u*u) + (B_v*v)
%
% Discrete time:
% x_k+1 = (Phi*x_k) + (Gamma_u*u_k) + (Gamma_v*v_k)
%
% This conversion is accomplished with the dynamics subfunction and ode
% solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) x - the (nx x 1) state vector at time k
% 2) u - the (nu x 1) disturbance input vector at time k.  It is assumed 
%       to be held constant throughout the integration interval.
% 3) v - the (nv x 1) control input vector at time k.  It is 
%       assumed to be held constant throughout the integratio interval.
% 4) tspan- the time span over which to evaluate the dynamics.  This
%       can be a two element vector [t_start, t_end] or a column vector.
%       The units are assumed to be in seconds.
% 5) param - this is a parameter structure that contains all of the other
%       inputs that might be needed by various dynamics models. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1) t - the final time that the state is defined at, in seconds.
%   2) x - the (nx x 1) state vector defined at time tspan(end). 
%   3) Phi - the (nx x nx) state transition matrix.  If it were a linear 
%           system, then the STM would take a state from time t_k and 
%           propagates it to time t_kp1.  It is just an approximation.
%   4) Gamma_u - the (nx x nu) process noise influence matrix that takes a
%           disturbance input vector and defines its influence from 
%           t_k->t_kp1.
%   5) Gamma_v - the (nx x nv) control influence matrix that takes a
%           control input vector and defines its influence from t_k->t_kp1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) The standard assumptions for ode solvers.  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x] = Provided_dynFxn_2Body_ECI_stateOnly(x0,tspan,param)

% ODE solver parameters:
odeSolver = @ode113;
dynarg    = param;
% ODE solver parameters:
dynarg    = param;
if isfield(param, 'odeSolver')
    odeSolver = param.odeSolver;
else
    odeSolver = @ode113;
end    
if isfield(param, 'odeOpts')
    odeOpts = param.odeOpts;
else
    %ODE Function Options
    odeOpts = odeset('abstol',3e-14,'reltol',3e-14,'MaxStep',100); %ODE Function Options
end
%
%% Form Initial State/STM
nx   = length(x0);
nx2  = nx^2;
z0   = zeros(nx,1);
% Populate State:
z0(1:nx,:) = x0;
%
%% Evaluate ODE:
odefunCall = @(t,x) odefun(t,x,@dynFxn_ContinuousTime_2Body_ECI,nx,dynarg);
[t,z]      = feval(odeSolver,odefunCall,tspan,z0,odeOpts);
%
%% Expand/Unflatten Solution
x       = z(:,1:nx)';
%
% Solver may return values at multiple times, but we only want the last
% one:
% x       = x(:,end);
end

function zdot = odefun(t,z,dynfun,nx,dynarg)
%% Sizing:
nx2  = nx^2;
%% Parse  inputs and unflatten (vector->matrix):
% Grab state:
x = z(1:nx,:);
%
%% Call Dynamics:
[xdot] = feval(dynfun,t,x,0,0,dynarg);
%
%% Flatten (matrix->vector):
% State:
zdot(1:nx,1)  = xdot;
end



%% Setup function to be called by the ode solver:
function [dx] = dynFxn_ContinuousTime_2Body_ECI(t,x,u,v,param)

% Earth mu:
muE_m = 3.986004415e14;
%
R = norm(x(1:3));
% Sanity check rho_km:
if (abs(R)<6371e3)
    error('Object collided with the Earth.')
end
%
% Distances in x,y,z:
RX = x(1,:);
RY = x(2,:);
RZ = x(3,:);
%
% Compute Gravitational Acceleration
%
% Define Constants
C0 = -muE_m./R.^3;
% Define Additional Constants
C2 = 3*muE_m./R.^5;
%
xdd = C0.*RX;
ydd = C0.*RY;
zdd = C0.*RZ;
dd  = [xdd;ydd;zdd];
%
% Form State Derivative 
% dx = zeros(size(x));
nStates   = size(x,1);
dx        = zeros(nStates,1);
dx(1:3,1) = x(4:6); %d/dt x = xdot 
dx(4:6)   = dd(:); % d/dt xdot = acceleration
end
