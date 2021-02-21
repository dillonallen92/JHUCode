clear, clc, close all;

% Constants/values we may need 

mu = 132712440041.94; % km^3/s^2
AUtoKm = 149597870.7; % km 
a = 1.32489106722386 * AUtoKm; % km
e = 0.255915584525353; % unitless
rp = 0.98583079532113; % AU
ra = 1.66395133912701; % AU
Period = 557.017301288787; % Days
p = a*(1-e^2); % km 

%% Problem 4

% Calculate the amount of time that starman spends its orbit under 1 AU.
% Note that you must make sure that the second mean anomaly is greater than
% the first by adding 2*pi

% Angular drift is n = 2pi/Period

n = 2*pi/Period; % radians / day
r_1AU = 1*AUtoKm; % km
theta_1AU = acos((p/r_1AU - 1)/e); % radians

E_1AU = 2*atan(sqrt((1-e)/(1+e))*tan(theta_1AU/2)); % radians 
M_1AU = E_1AU - e*sin(E_1AU); % radians 
t_1AU = M_1AU / n; % days

% Because of symmetry, we can multiply this value by 2 to get the full time
% under 1 AU

t_1AU_total = 2 * t_1AU; % days


%% Problem 5

% Independently calculate the amount of time that starman spends above 1
% AU and verify that this and the value obtained in 4 add up to the orbital
% period.

% The approach I will do for this problem will be to calculate the time it
% takes to get to 1AU on the first side of the ellipse and then do 2*pi -
% theta_1AU to find M_prime and t_prime at the other side. The difference
% in time should be the total time spent above 1AU through the elliptical
% orbit.

% Using Q4 as an abbreviation for Quadrant 4

theta_1AU_Q4 = 2*pi - theta_1AU; % radians 
E_1AU_Q4 = 2*atan(sqrt((1-e)/(1+e))*tan(theta_1AU_Q4/2)); % radians 

% The value should come out negative because we are in Q4 and tan gives us
% the value of the reference angle, so add 2*pi to get the angle in Q4

if E_1AU_Q4 < 0
    E_1AU_Q4 = E_1AU_Q4 + 2*pi; % radians
end

M_1AU_Q4 = E_1AU_Q4 - e*sin(E_1AU_Q4); % radians 
t_1AU_Q4 = M_1AU_Q4 / n; % days

% Time spent outside of 1 AU
t_outside1AU = t_1AU_Q4 - t_1AU; % days

% Confirm they add up to 557.0173 days
t_tot = t_1AU_total + t_outside1AU; % days
