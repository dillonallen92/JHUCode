clear, clc, close all;

% Weighted Least Squares Example as seen in slide 31
% Dynamics model provided: r = x1 + x2 * t;
% Measurements and noise distribution provided
% Nosie = w ~ N(0,0.2*i)

y = [ .99 1.93 3.37 4.87 6.1]';
t = 1:5;
H = [1 0; 1 1; 1 2; 1 3; 1 4];
v = (0.2*t).^2;
R = diag(v);

% No Prior info
x_hat = inv(H' * inv(R) * H)*H'*inv(R)*y;
P = inv(H' * inv(R) * H);

% Compare to Least Squares
x_hat_LS = inv(H'*H)*H'*y;