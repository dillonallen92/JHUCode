clear, clc, close all;

% Use the following state vector 
xvec = [7378.05e3; -33.75; -14.667; 0.03; 7.35e3; 0];

%% Part A
% for zero control and disturbance, prop for 105 min.
tspan = 0:0.1:(105*60);
u = [ 0 0 0]';
v = [ 0 0 0]';
param = 0;
[t,x,Phi,Gamma_u,Gamma_v] = dynFxn_2Body_ECI_Template(xvec,u,v,tspan, param);
