clear, clc, close all;

% Problem 3
% You are designing a sun synchronous orbit 
% (Omega_dot = 2pi/365.25 rad/day) with an eccentricity of 0.3 and a
% semimajor axis of 10,000 km

a = 10000; % km
e = 0.3;
Omega_dot = 2*pi/365.25 * (1/(24*3600)); % rad/sec
r = 6378; % km
J2 = 0.001082;
mu_E = 3.986e5; % km^3/s^2

n = sqrt(mu_E/a^3); % rad / s
p = a*(1-e^2); % km

%% Part A
% Find the inclination i

i = acos(-2*Omega_dot*p^2 / (3*n*J2*r^2));

%% Part B
% find the average change of the argument of periapsis
omega_dot = (3/4)*n*J2*(r/p)^2 * (5*cos(i)^2 - 1);