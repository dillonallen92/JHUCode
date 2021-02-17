clear, clc, close all;

% Problem 6
% Calculate the delta-v necessary to put starman on an escape trajectory
% with respect to the sun at any point along his orbit. Make a plot of
% delta-v vs true anomaly. Where is the minimum value? Why does this make
% sense?

% delta-v formula: dv(r(f)) = sqrt(2*mu / r(f)) - v(r(f))
% v(r(f)) = sqrt( mu * (2/r(f) - 1/a))
% r(f) = (a*(1-e^2))/(1+e*cos(f))

% Constants (a and e from problem 1)
mu = 132712440041.94 % km^3/s^2
AUtoKm = 149597870.7 % km 
a = 1.32489106722386 * AUtoKm; % km
e = 0.255915584525353; % unitless

f = 0:0.01:pi; % rads
r = (a*(1-e^2))./(1+e*cos(f)); % km
v = sqrt( mu * (2./r - 1/a)); % km/s
dv = sqrt(2*mu ./ r) - v ; % km/s
plot(f, dv);
title("\DeltaV vs Mean Anomaly");
xlabel("Mean Anomaly (rad)");
ylabel("\DeltaV (km/s)");