clear, clc, close all;

% Find the time of flight between r1 and r3

% Position vectors about the sun
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

deltaT = sqrt(a^3/mu_S)*(alpha - beta - (sin(alpha) - sin(beta)));
days = deltaT/(24*3600);