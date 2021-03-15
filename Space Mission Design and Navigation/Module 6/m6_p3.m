clear, clc, close all;

% If instead of a Hohmann transfer we decide to put the sc on a transfer orbit with
% radius of periapsis equal to the radius of Venus' orbit and radius of apoapsis equal
% to 10 times the radius of Saturns orbit.

% Values that may be important
rE = [141725773.204927; -51913278.0968566; 2454.00560975075]; % km
vE = [9.75958833607707; 27.8458269687754; -0.00137645277890108]; % km/s
rS = [-690232102.803463; 1168541543.04847; 7181595.50405902]; % km
vS = [-8.84280567480198; -4.93025746647439; 0.43850160809551]; % km/s

mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

rEnorm = norm(rE); % km
rSnorm = norm(rS); % km
rVnorm = 108.209e6; % km

%% Part A
% Calculate the tof, dv1, dv2

rp = rVnorm;
ra = 10*rSnorm;

a = (rp + ra)/2;

tof = pi*sqrt(a^3/mu_S);

dv1 = sqrt(mu_S*(2/rp - 1/a)) - sqrt(mu_S/rp);

dv2 = sqrt(mu_S/ra) - sqrt(mu_S*(2/ra - 1/a));

%% Part B
% What is the estimate for launch C3?

En = dv1^2/2;
C3 = 2*En;

%% Part C
% What is the new estimate for a dv required to capture into a 200 day orbit about Saturn?
mu_Saturn =3.7931187E16 *(1/1000)^3 ; % km^3/s^2
rp_Saturn = 63281.4; % km
P = 200 * 24 * 3600; % s
a_saturnOrbit = (mu_Saturn*P^2/(2*pi)^2)^(1/3); % km
E_Saturn = -mu_Saturn/(2*a_saturnOrbit);

v_SaturnOrbit = sqrt(2*(E_Saturn + mu_Saturn/rp_Saturn));
v_SC_Saturn = sqrt(2*(En + mu_Saturn/rp_Saturn));
dv_SaturnOrbit = v_SaturnOrbit - v_SC_Saturn; % 