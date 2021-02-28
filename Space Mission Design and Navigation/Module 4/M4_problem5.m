clear, clc, close all;

% A 500kg spacecraft in a circular 350km altitude orbit will need to
% perform an orbit maintenance maneuver each time that its altitude drops
% below 340 km to complete its mission. Assume that Cd = 2.0, A = 2 m^2 and
% the density is constant over this range, rho = 1.225E-11 kg/m^3. How long
% will the spacecraft take to drop 10km before its stationkeeping maneuver?

% Problem Constants
m = 500; % kg
Ri = 350; % km
Rf = 340; % km
Cd = 2.0;
A = 2; % m^2
rho = 1.225E-11; % kg/m^3
mu = 3.986004418E14 * (1/1000)^3 ; % km^3/s^2

Bstar = Cd * A / m;

time = -2 * (Rf^(1/2) - Ri^(1/2)) / (mu^(1/2) * rho * Bstar); % seconds
days = time / (3600 * 24);