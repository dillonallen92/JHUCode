clear, clc, close all;

% Planar Two-Body problem
% x_double_dot = -mu X/r^3
% y_double_dot = -mu Y/r^3
% xvec = [X Y Vx Vy]

%% Part A
% Linearize the system
% r = sqrt(xvec(1)^2 + xvec(2)^2
% xvec_dot = [xvec(3) xvec(4) -mu*xvec(1)/r^3 -mu*x(2)/r^3];

%% Part B
% Derive each element of the A Matrix
% A = [ 0                        0                        1 0; 
%       0                        0                        0 1;
%       3*x(1)^2*mu/r^5 - mu/r^3 3*x(1)*x(2)*mu/r^5       0 0;
%       3*x(1)*x(2)*mu/r^5       3*x(2)^2*mu/r^5 - mu/r^3 0 0;]   
%

%% Part C
% Write a matlab function for integration which combines the equations of
% motion. The function will have the form
% [xdot, t] = twoBody(t,x,mu)
% and called via
% [t, x] = ode113(@(t,x)twoBodySTM(t,x,mu),tspan,xo,options)
% x is a 20x1 vector
% Use the reshape command,
% phi = transpose(reshape(x(5:20,1),4,4));
% phi_dot = A * phi;
% x_dot(5:20,1) = reshape(transpose(phi_dot),16,1);

%% Part D
% The spacecraft is in a circular orbit with the following initial state
% about earth
r = [42241.0800678832 0]'; % km
v = [0 3.07185802826297]'; % km/s
muE = 3.986e5;
tspan = 0:10:12*3600;
options=odeset('abstol',1e-10,'reltol',1e-10);
xComb0 = [r(1); r(2); v(1); v(2); reshape(eye(4),16,1)];
[t,x] = ode113(@(t,x)twoBody(t,x,muE),tspan,xComb0, options);

xfinal = x(end,1:20);

%% Part E
% Reform the last state back into the 4x4 STM
% Identify the matrix element which predicts the change in the final
% x-component of position with the initial y-component of velocity. Use
% this element to predict the change in the final x-component of position
% for a 0.001 km/s increase in Vy.

STMFinal = reshape(xfinal(1,5:20),4,4)';

%% Part F
% Propogate with part d
rnew = [42241.0800678832 0]';
vnew = [0 3.07285802826297]';
xCombNew = [rnew(1); rnew(2); vnew(1); vnew(2); reshape(eye(4),16,1)];
[tnew, xnew] = ode113(@(t,x)twoBody(t,x,muE),tspan,xCombNew,options);

xnew_final = xnew(end, 1:20);
xnew_final_x = xnew(1,1);
xfinal_x = xfinal(1,1);
partF_difference = xnew_final_x - xfinal_x;


%% Functions
% Part C Function

function [x_dot, t] = twoBody(t,x,mu)
    r = sqrt(x(1)^2 + x(2)^2);
    x_dot(1:4,1) = [x(3) x(4) -mu*x(1)/r^3 -mu*x(2)/r^3]';
    A = [ 0                        0                        1 0; 
          0                        0                        0 1;
          3*x(1)^2*mu/r^5 - mu/r^3 3*x(1)*x(2)*mu/r^5       0 0;
          3*x(1)*x(2)*mu/r^5       3*x(2)^2*mu/r^5 - mu/r^3 0 0;];
    phi = transpose(reshape(x(5:20,1),4,4));
    phi_dot = A*phi;
    x_dot(5:20,1) = reshape(transpose(phi_dot),16,1);
end

