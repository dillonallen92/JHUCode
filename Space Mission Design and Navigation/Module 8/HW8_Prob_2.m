clear, clc, close all;

%% Part 1
% Load the input from HW8_Prob2.mat. You will need uHist and dt.
load('HW8_Prob2.mat');

%% Part 2
% Simulate the clock drift for 1000 cycles of dt, using the new process
% noise at each step. Plot the state time-history output

tspan = 0:dt:1000*dt;
v = [0 0 0]';
x = [0 0 0]';
u = uHist;
param.derFlag = 1;
[t,x,Phi,Gamma_u,Gamma_v] = dynamicsFunc(x,u,v,tspan,param);

% Since I used a function handle, I plugged in the last two time steps to
% calculate the state transition matrix at the end
phi_mat = Phi(tspan(length(tspan)), tspan(length(tspan)-1));

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);


subplot(311)
plot(tspan,x1);
title('Time History Plots');
xlabel('time (s)');
ylabel('Clock Error');
subplot(312)
plot(tspan,x2);
xlabel('time (s)');
ylabel('Frequency Error');
subplot(313)
plot(tspan,x3);
xlabel('time (s)');
ylabel('Frequency Rate Error');