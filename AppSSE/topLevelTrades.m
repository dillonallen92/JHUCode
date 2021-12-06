% top level trades to Jupiter
% Delta-V calculation and C3/mass ratio

clear, clc, close all;

% Constants
m0 = 588;   % kg
Isp = 315;  % s
dV = 705;   % m/s
T = 150;    % lbf
g0 = 9.81;   % m/s^2

% Convert lbf to N
T_N = 150 * 4.44822;    % Newtons

% Exhaust Velocity
ve = g0*Isp;            % m/s

% Final mass
mf = m0*exp(-dV/ve);    % kg

dt = (-ve*(mf-m0))/T_N; % s
dt_min = dt/60;