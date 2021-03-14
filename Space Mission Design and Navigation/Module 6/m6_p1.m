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