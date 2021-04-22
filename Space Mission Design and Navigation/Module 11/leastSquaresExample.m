clear, clc, close all;

% Least Squares Example as seen in slide 14
% Dynamics model: r = x1 + x2 * t;
% Measurement State Provided
% t = 0 - 4

% Measurement State
y = [ 1.1 2.36 2.54 4.17 5.06]';

% Jacobian
% H = [ 1 t]
H = [ 1 0; 1 1; 1 2; 1 3; 1 4];

% estimated state
x_hat = inv(H' * H) * H' * y;
x_true = [ 1 1 ]';