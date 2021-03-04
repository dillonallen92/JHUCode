clear, clc, close all;

% Solve the transfer with Lambert's Problem

% Info from problem 3
r1 = [-131386230.977293, 69971484.9501445, -718889.822774674]'; % km
r3 = [-139952788.024352, 54396927.5494364, -1043119.92803938]'; % km

% Sun gravitational parameter
mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

% Ask the professor but I believe I can use my a from problem 2

r1norm = norm(r1);
r3norm = norm(r3);
a = 1.982008825664394e+08;
c = norm(r3-r1);
s = (r1norm + r3norm + c)/2;
alpha = 2*asin(sqrt(s/(2*a)));
beta = 2*asin(sqrt((s-c)/(2*a)));

%% Part A
% Find the minimum semi-major axis and the corresponding transfer time to
% go from r1 to r3

a_min = s/2;
beta_min = 2*asin(sqrt((s-c)/(s)));

tmin = sqrt(s^3 / (8*mu_S))*(pi - beta_min + sin(beta_min));
minDays = tmin*(1/(3600*24));
