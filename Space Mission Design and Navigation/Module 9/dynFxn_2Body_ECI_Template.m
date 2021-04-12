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
% Difference Equation:
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
function [t,x,Phi,Gamma_u,Gamma_v] = dynFxn_2Body_ECI_Template(x0,u,...
    v,tspan,param)

% ODE solver parameters:
odeSolver = @ode113;
dynarg    = param;
odeOpts   = odeset('abstol',1e-13,'reltol',1e-13,'MaxStep',1000); % Recommended.
%
%% Form Initial State/STM
nx   = length(x0);
nx2  = nx^2;
nu   = length(u);
nv   = length(v);
nvnx = (nx*nv);
nunx = (nx*nu);
z0   = zeros((nx+nx2+nunx+nvnx),1);
% Populate State:
z0(1:nx,:) = x0;
% Populate STM:
initial_dxdotdx   = eye(nx);
z0(nx+1:nx+nx2,:) = reshape(initial_dxdotdx,nx2,1);
% Populate CIM:
initial_dxdotdu   = zeros(nx,nu);
z0((nx+nx2)+1:(nx+nx2+nunx),:) = reshape(initial_dxdotdu,nunx,1);
% Populate PNIM:
initial_dxdotdv   = zeros(nx,nv);
z0((nx+nx2+nunx)+1:(nx+nx2+nunx)+nvnx,:) = reshape(initial_dxdotdv,nvnx,1);
%
%% Evaluate ODE:
odefunCall = @(t,x) odefun(t,x,@dynFxn_ContinuousTime_2Body_ECI,nx,nu,u,nv,v,dynarg);
[t,z]      = feval(odeSolver,odefunCall,tspan,z0,odeOpts);
%
%% Expand/Unflatten Solution
nt      = length(t);
x       = z(:,1:nx)';
Phi     = reshape(z(:,nx+1:nx+nx2)',nx,nx,nt);
Gamma_u = reshape(z(:,(nx+nx2+1):(nx+nx2+nunx))',nx,nu,nt);
Gamma_v = reshape(z(:,(nx+nx2+nunx+1):end)',nx,nv,nt);
%
% Solver may return values at multiple times, but we only want the last
% one:
x       = x(:,end);
Phi     = Phi(:,:,end);
Gamma_u = Gamma_u(:,:,end);
Gamma_v = Gamma_v(:,:,end);
end

function zdot = odefun(t,z,dynfun,nx,nu,u,nv,v,dynarg)
%% Sizing:
nx2  = nx^2;
nunx = nu*nx;
nvnx = nv*nx;
%% Parse  inputs and unflatten (vector->matrix):
% Grab state:
x = z(1:nx,:);
% Grab STM:
Phi     = reshape(z(nx+(1:nx2)),nx,nx);
% Expand (unflatten) Process Noise Influence Matrix:
Gamma_u = reshape(z((nx+nx2)+1:(nx+nx2+nunx)),nx,nu);
% Expand (unflatten) Control Influence Matrix:
Gamma_v = reshape(z((nx+nx2+nunx)+1:(nx+nx2+nunx+nvnx)),nx,nv);
%
%% Call Dynamics:
[xdot,A,Bu,Bv] = feval(dynfun,t,x,u,v,dynarg);
%
%% Flatten (matrix->vector):
% State:
zdot(1:nx,1)                 = xdot;
% Phi:
Phidot                       = A * Phi;
zdot(nx+(1:nx2),1)           = Phidot(:);
% Control and Process Noise:
Gamma_uDot                   = A*Gamma_u + Bu;
zdot(nx+nx2+(1:nunx),1)      = Gamma_uDot(:);
Gamma_vDot                   = A*Gamma_v + Bv;
zdot(nx+nx2+nunx+(1:nvnx),1) = Gamma_vDot(:);
end



%% Setup function to be called by the ode solver:
function [dx,A,Bu,Bv] = dynFxn_ContinuousTime_2Body_ECI(t,x,u,v,param)

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
% Process Noise:
uX = u(1,:);
uY = u(2,:);
uZ = u(3,:);
%
% Control:
vX = v(1,:);
vY = v(2,:);
vZ = v(3,:);
%
% Compute Gravitational Acceleration
%
% Define Constants
C0 = -muE_m./R.^3;
% Define Additional Constants
C2 = 3*muE_m./R.^5;
%
xdd = C0.*RX + vX + uX;
ydd = C0.*RY + vY + uY;
zdd = C0.*RZ + vZ + uZ;
dd  = [xdd;ydd;zdd];
%
% Form State Derivative 
% dx = zeros(size(x));
nStates   = size(x,1);
dx        = zeros(nStates,1);
dx(1:3,1) = x(4:6); %d/dt x = xdot 
dx(4:6)   = dd(:); % d/dt xdot = acceleration
%
%% Construct Jacobian (if called)
% Partial Derivative Matrix w.r.t. state:
A(1:3,4:6) = eye(3); %dx1dx1
%
A(4,1)     = C0 + C2.*RX.^2; %dx2dx
A(4,2)     =      C2.*RX.*RY;%dx2dy
A(4,3)     =      C2.*RX.*RZ;%dx2dz
%
A(5,1)     =      C2.*RY.*RX;%dy2dx
A(5,2)     = C0 + C2.*RY.^2; %dy2dy
A(5,3)     =      C2.*RY.*RZ;%dy2dz
%
A(6,1)     =      C2.*RZ.*RX;%dz2dx
A(6,2)     =      C2.*RZ.*RY;%dz2dy
A(6,3)     = C0 + C2.*RZ.^2; %dz2dz
%
% Partial Derivative Matrix w.r.t. control noise:
Bu = zeros(nStates,3);
Bu(4,1) = 1; %dx2dux
Bu(5,2) = 1; %dy2duy
Bu(6,3) = 1; %dz2duz
%
% Partial Derivative Matrix w.r.t. process noise:
Bv = zeros(nStates,3);
Bv(4,1) = 1; %dx2dvx
Bv(5,2) = 1; %dy2dvy
Bv(6,3) = 1; %dz2dvz
end
