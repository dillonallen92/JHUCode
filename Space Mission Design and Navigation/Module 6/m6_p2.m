clear, clc, close all;

% Assume that Earth and Saturn are on circular orbits. Neglecting their gravity, calculate
% the Hohmann transfer to go from Earth to Saturn

% Values from previous problem
rE = [141725773.204927; -51913278.0968566; 2454.00560975075]; % km
vE = [9.75958833607707; 27.8458269687754; -0.00137645277890108]; % km/s
rS = [-690232102.803463; 1168541543.04847; 7181595.50405902]; % km
vS = [-8.84280567480198; -4.93025746647439; 0.43850160809551]; % km/s

mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

rEnorm = norm(rE); % km
rSnorm = norm(rS); % km
tof_days = 1214; % days

%% Part A
% Calculate the time of flight, dv1 and dv2

a_h = (rEnorm + rSnorm) / 2; % km
tof = pi*sqrt(a_h^3/mu_S); % s

E_h = -mu_S/(2*a_h); % kJ

v_pH = sqrt(2*(E_h + mu_S/rEnorm)); % km / s
v_c1 = sqrt(mu_S/rEnorm); % km /s

dv1 = v_pH - v_c1; % km/s 

v_c2 = sqrt(mu_S/rSnorm); % km/s
v_aH = sqrt(2*(E_h + mu_S/rSnorm)); % km/s

dv2 = v_c2 - v_aH; % km / s

%% Part B
% A launch from Earth can be estimated assuming v_inf = dv1. What is the C3 and how does
% this value ompare to the C3 of the truncated trajectory from the previous problem?

C3_p1 = 73.5543;
C3 = dv1^2;

dC3 = C3_p1 - C3;

%% Part C
% What would be our new estimate for a dv required to capture into a 200 day orbit about 
% Saturn if we make the approximation that v_inf = dv2 for Saturn?

v_inf = dv2;
En = v_inf^2/2;

mu_Saturn =3.7931187E16 *(1/1000)^3 ; % km^3/s^2
rp_Saturn = 63281.4; % km
P = 200 * 24 * 3600; % s
a_saturnOrbit = (mu_Saturn*P^2/(2*pi)^2)^(1/3); % km
E_Saturn = -mu_Saturn/(2*a_saturnOrbit);

v_SaturnOrbit = sqrt(2*(E_Saturn + mu_Saturn/rp_Saturn));
v_SC_Saturn = sqrt(2*(En + mu_Saturn/rp_Saturn));
dv_SaturnOrbit = v_SaturnOrbit - v_SC_Saturn; 

%% Part D
% Given the tof what is the relative angular position between Earth and Saturn at launch?

P_Saturn = (2*pi/(sqrt(mu_Saturn)))*(rSnorm)^(3/2);
halfP_days = P_Saturn/2 * (1/(24*3600));
deg_per_day = 180/halfP_days;
angular_sep = (halfP_days - tof_days) * deg_per_day;
