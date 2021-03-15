clear, clc, close all;

m0 = 1200; % kg
g0 = 9.81; % m/s^2
Isp = 225; % s
thrust = 22; % N

% Maneuvers (m/s)
l_cleanup = 20;
v1_targeting = 10;
v1_cleanup = 5;
e1_targeting = 10;
e1_cleanup = 5;
e2_targeting = 10;
s1_targeting = 20;

%% Part A
% If the maneuvers include margin, how much propellant in kg should be allocated

maneuver_list = [l_cleanup, v1_targeting, v1_cleanup, e1_targeting, e1_cleanup, ...
                           e2_targeting, s1_targeting];

prop_mass = 0;
dm = 0;
mi = m0;
for i = 1 : length(maneuver_list)
    m1 = mi - dm;
    dv = maneuver_list(i);
    dm = m1*(1-exp(-(dv)/(Isp*g0)));
    prop_mass = prop_mass + dm;
end

% testing function
prop_mass1 = maneuver_prop_calc(m0, maneuver_list, g0, Isp);

fprintf("Total propellant mass: %d\n", prop_mass1);

%% Part B
% Make a plot that shows how the propellant increases with mass to 5,000 kg

m_final = 5000; % kg
mList = m0:1:5000; % kg
propList = maneuver_prop_calc(mList, maneuver_list, g0, Isp);

% Plot of Mass vs Propellant
plot(mList, propList);
xlabel(" SC Mass (kg) ");
ylabel(" Propellant Allocation (kg) ");
title(" SC mass vs Propellant Allocation ");

%% Part C
% What is the mass flow rate for this engine?
mdot = thrust/(Isp * g0);

%% Part D
% The positions and velocity of Earth at E2 and Saturn at S1 in the ecliptic J2000 frame are
% below. Reconstruct the final leg of this trajectory with the lambert solver given tof = 1214 days.
% List the value of semi-major axis, initial and final velocities.

rE = [141725773.204927; -51913278.0968566; 2454.00560975075]; % km
vE = [9.75958833607707; 27.8458269687754; -0.00137645277890108]; % km/s
rS = [-690232102.803463; 1168541543.04847; 7181595.50405902]; % km
vS = [-8.84280567480198; -4.93025746647439; 0.43850160809551]; % km/s

mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

rEnorm = norm(rE); % km
rSnorm = norm(rS); % km
c = norm(rS-rE); % km
s = (rEnorm + rSnorm + c)/2; % km

a_min = s/2; % km 
beta_min = 2*asin(sqrt((s-c)/(s))); % radians

tmin = sqrt(s^3 / (8*mu_S))*(pi - beta_min + sin(beta_min)); % s
minDays = tmin*(1/(3600*24)); % days

tof = 1214 * 24 * 60 * 60; % seconds

theta = acos(dot(rE,rS)/(rEnorm*rSnorm)); % radians

[a_sol, fval, exitflag, output] = fzero(@(x) lamFun(x,tof,s,c,theta,tmin,mu_S), [a_min, 10*a_min]);

alpha_vel = 2*asin(sqrt(s/(2*a_sol)));
beta_vel = 2*asin(sqrt((s-c)/(2*a_sol)));
A = sqrt(mu_S/(4*a_sol))*cot(alpha_vel/2); % km/s
B = sqrt(mu_S/(4*a_sol))*cot(beta_vel/2); % km/s 
uc = (rS-rE)/norm(rS-rE); 
uE = rE/rEnorm;
uS = rS/rSnorm;

vi = (B+A)*uc + (B-A)*uE; % km/s

vf = (B+A)*uc - (B-A)*uS;

fprintf("Semi-major axis of trajectory: %d (km)\n", a_sol);
fprintf("Initial Velocity: ");
disp(vi);
fprintf("Final Velocity: ");
disp(vf);

%% Part E
% Given the final velocity of the spacecraft w.r.t the sun, estimate the excess velocity
% w.r.t Saturn and calculate the energy w.r.t Saturn

v_inf = vf - vS;

En = dot(v_inf, v_inf) / 2;

%% Part F
% If a periapsis radius of 63,281.4 km is targeted, calculate the dv required to capture into
% orbit with period of 200 days

mu_Saturn =3.7931187E16 *(1/1000)^3 ; % km^3/s^2
rp_Saturn = 63281.4; % km
P = 200 * 24 * 3600; % s
a_saturnOrbit = (mu_Saturn*P^2/(2*pi)^2)^(1/3); % km
E_Saturn = -mu_Saturn/(2*a_saturnOrbit);

v_SaturnOrbit = sqrt(2*(E_Saturn + mu_Saturn/rp_Saturn));
v_SC_Saturn = sqrt(2*(En + mu_Saturn/rp_Saturn));
dv_SaturnOrbit = v_SaturnOrbit - v_SC_Saturn; % should be negative to indicate slowing down

m0_atSaturn = m0 - prop_mass;
saturn_prop = maneuver_prop_calc(m0_atSaturn, abs(dv_SaturnOrbit), g0, Isp);

time_for_maneuver = saturn_prop / mdot; % seconds

m_remaining = m0_atSaturn - saturn_prop;

%% Part G
% C3 = 2*En. If C3 ~ 16 km^2/s^2, and they cancel the inner legs, what is the difference
% in C3?

% C3 from E2 to S1
C3_E2S1 = 2*En;
C3_problemStatement = 16;
dC3 = C3_E2S1 - C3_problemStatement;

% Discussion: Removing the inner legs of the trip will require more direct launch energy
% to accomplish the mission, since we cannot rely on gravity assists. This will make the
% mission harder in the fact that you have to find rockets able to match the higher C3 
% needed as well as match the needed SC mass.


%% Functions

function totalProp = maneuver_prop_calc(m0, maneuver_list, g0, Isp)
    prop_mass = 0;
    dm = 0;
    mi = m0;
    for i = 1 : length(maneuver_list)
        m1 = mi - dm;
        dv = maneuver_list(i);
        dm = m1*(1-exp(-(dv)/(Isp*g0)));
        prop_mass = prop_mass + dm;
    end
    
    totalProp = prop_mass;
end

function val = lamFun(x,tof,s,c,theta,tmin, mu_S)

    alpha_0 = 2*asin(sqrt(s/(2*x)));
    beta_0 = 2*asin(sqrt((s-c)/(2*x)));

    if (0 <= theta) && (theta <= pi)
        beta = beta_0;
    else
        beta = - beta_0;
    end
    
    if tof <= tmin
        alpha = alpha_0;
    else
        alpha = 2*pi - alpha_0;
    end
    
    val = sqrt(mu_S/x^3)*tof -  (alpha - beta - (sin(alpha) - sin(beta)));
end
