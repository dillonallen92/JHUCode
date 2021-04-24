clear, clc, close all;

% Attempt on the Az/El example as seen on slides 37 - 41
load('NLWLS_Example.mat');

% Process
% Convert from ECEF to ENV Coordinates
% Compute az/el
% Find Jacobian through numerical differencing
% x*_0 = x_0
% x_hat(1) = x_0 - inv(H'(x_0)*H(x_0))*(H'(x_0)*(h(x_0) - z(1)))
% x_star_i+1 = x_hat(i)
% repeat until norm(x_star(i+1) - x_star(i)) < 1e-6

% Questions:
% where do I get the measurements for the example problem
% How do I do this problem with 3 instrument locations
% The flow chart in 4.6 uses the STM, but the notes have a process without it, which one?