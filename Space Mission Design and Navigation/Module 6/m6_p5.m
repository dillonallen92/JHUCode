clear, clc, close all;

% A sc is maneuvering from a circular LEO of 200 km to a circular orbit of a much greater
% altitude of 92294.1 km
rE = 6378; % km
rp = 200 + rE;
ra = 92294.1 + rE;

mu_E = 3.986004418E14 * (1/1000)^3; % km^3/s^2

%% Part A
% Calculate the tof, dv1, dv2 for the Hohmann Transfer

a_H = (rp + ra)/2;

tof = pi*sqrt(a_H^3/mu_E);

dv1 = sqrt(mu_E*(2/rp - 1/a_H)) - sqrt(mu_E/rp);

dv2 = sqrt(mu_E/ra) - sqrt(mu_E*(2/ra - 1/a_H));

%% Part B
% Calculate dv1, dv2, dv3 and the tof for a bi-elliptic transfer for ri = 2r2.

ri = 2*ra;
a_new = (ri + rp)/2;



%% Part C
% As discussed in the textbook the total dv for the bi-elliptic transfer should be less than
% the Hohmann option. What are some of the limiting factors that could forbid Earth orbiting 
% sc from executing a bi-elliptic transfer for even higher altitudes