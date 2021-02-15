clear, clc, close all;

% The system above will be used in this problem but in the 
% J2000 Frame, as a single mass around the sun

% Given variables and constants
r = [159534320.691; -57419276.2529436; -10082710.3505116];
v = [3.33434967532767; 31.8108643677972; 0.357944743421687];
mu = 1.327e11;


%% Part A
% Calculate its period about the sun

% First I need to find the semi-major axis, so I will use a 
% technique I did from problem 1

E = dot(v,v)/2 - mu/norm(r);
a = - mu / (2 * E);

Period = 2*pi*sqrt(a^3/mu);
Period_days = Period / (24 * 3600);

%% Part B
% Calculate the eccentricity 

% Also pulling from the equation I used in problem 1 to find
% the eccentricity

h = cross(r,v);
h_norm = norm(cross(r,v));
e = sqrt(1 - h_norm^2/(mu * a));

%% Part C
% What is the radius of periapsis and apoapsis?

rp = a*(1-e);
ra = a*(1+e);

%% Part D
% What is the angular momentum vector?

% Well, I already calculated this above to find my eccentricity

%% Part E
% What is the time rate of change of the truye anomaly and
% the time rate of change of the Area swept out by its radius
% vector with respect to the sun?

% the rate of change of the true anomaly will be from 
% h = r^2 * theta_dot

r_norm = norm(r);
theta_dot = h_norm / r_norm^2;

dA_dt = h_norm / 2;

%% Part F
% Calculate its position and velocity in the perifocal frame.

% First I need to calculate the mean anomaly with the given
% r, a and e.

p = h_norm^2 / mu;
theta = acos((p/r_norm - 1)/e);
v_norm = norm(v);

r_peri = r_norm* [cos(theta); sin(theta); 0];
v_peri = [v_norm * cos(theta) - r_norm * theta_dot * sin(theta);
            v_norm*sin(theta) + r_norm*theta_dot*cos(theta);
            0];

