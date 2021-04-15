clear, clc, close all;

% Use the following state vector 
xvec = [7378.05e3; -33.75; -14.667; 0.03; 7.35e3; 0]; % m ; m/s

%% Part 2
% for zero control and disturbance, prop for 105 min.
tspan = 0:0.1:(105*60);
u = [ 0 0 0]';
v = [ 0 0 0]';
param = 0;
[t,x,Phi,Gamma_u,Gamma_v] = dynFxn_2Body_ECI_Template(xvec,u,v,tspan, param);

%% Part 3
% Propagate the state forward for two cases for only 10 min 
tspan_new = 0:0.1:(10*60);
%  Case a
ua = [ 0 0 0]'; va = [0 0 0]';
% Case b
ub = [0.1 0.1 0.1]'; vb = [0.1 0.1 0.1]';

[ta, xa, Phia, Gamma_ua, Gamma_va] = dynFxn_2Body_ECI_Template(xvec, ua, va, tspan_new, param);
[tb, xb, Phib, Gamma_ub, Gamma_vb] = dynFxn_2Body_ECI_Template(xvec, ub, vb, tspan_new, param);

%% Part 4
% Provide the difference
diff = xa - xb;

%% Part 5
% My assumption here is f() is the dynFxn, and once we calculate that, we
% find the actual state x_a by f() + Phi_a * (xvec) + 0 (for case a)
[t_a, x_a, Phi_a, Gamma_u_a, Gamma_v_a] = dynFxn_2Body_ECI_Template(xvec, ua, va, tspan_new, param);
x_a = x_a + Phi_a * xvec + Gamma_u_a * [0 0 0]' + Gamma_v_a * [0 0 0]';

[t_b, x_b, Phi_b, Gamma_u_b, Gamma_v_b] = dynFxn_2Body_ECI_Template(xvec, ub, vb, tspan_new, param);
x_b = x_b + Phi_b*xvec + Gamma_u_a*ub + Gamma_v_a*vb;

approx_diff = x_a - x_b;

% There is a difference of like a factor of two so its either the linear
% approximation is "decent" or, the more likely case, I made some mistake
% and I cannot figure out what I did wrong. 