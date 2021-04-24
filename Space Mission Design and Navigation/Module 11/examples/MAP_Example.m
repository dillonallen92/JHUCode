clear, clc, close all;

% MAP Example
% Dynamics still the same: r = x1 + x2 * t
% same y and noise/covariance

y = [ .99 1.93 3.37 4.87 6.1]';
t = 1:5;
H = [1 0; 1 1; 1 2; 1 3; 1 4];
v = (0.2*t).^2;
R = diag(v);

% Prior Info
x_bar = [0.97; 1.1];
P_bar = [0.25^2 0; 0 0.25^2];

x_hat_MAP = inv(H' * inv(R) * H + inv(P_bar)) * (H' * inv(R) * y + inv(P_bar)*x_bar);
P_MAP = inv(H'*inv(R)*H + inv(P_bar));

% WLS
% No Prior info
x_hat_WLS = inv(H' * inv(R) * H)*H'*inv(R)*y;
P_WLS = inv(H' * inv(R) * H);

% Compare to Least Squares
x_hat_LS = inv(H'*H)*H'*y;