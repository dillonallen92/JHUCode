clear, clc, close all;

% Find the time of flight between r1 and r3

% Position vectors about the sun
r1 = [-131386230.977293, 69971484.9501445, -718889.822774674]'; % km
r3 = [-139952788.024352, 54396927.5494364, -1043119.92803938]'; % km

% Sun gravitational parameter
mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

r1norm = norm(r1); % km
r3norm = norm(r3); % km
a = 1.982008825664394e+08; % km
c = norm(r3-r1); % km
s = (r1norm + r3norm + c)/2; % km
alpha = 2*asin(sqrt(s/(2*a))); % radians
beta = 2*asin(sqrt((s-c)/(2*a))); % radians 

deltaT = sqrt(a^3/mu_S)*(alpha - beta - (sin(alpha) - sin(beta))); % s
days = deltaT/(24*3600); % days
