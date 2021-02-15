clear, clc, close all;

% Problem 1
% Suppose two masses are in free space with the given vector components

% Constants
G = 6.67E-17;

% Mass 1
m1 = 523026872797;
r1 = [10; 10; 0];
r1_dot = [0;1;0];

% Mass 2
m2 = 4773127203;
r2 = [11.1092062504984; 10; 0];
r2_dot = [0; 1.0001834740183; 0];

%% Part A
% Find the center of mass for the system

Rcm = (m1*r1 + m2*r2)/(m1+m2);

%% Part B
% What is the magnitude of the gravitational force in Newtons?
% Also what is it in lbs?
% Note to myself: The negative is just to denote attraction

F12 = -G*m1*m2/(norm(r1-r2))^2;

NewtonToLbConversion = 0.2248089431;

F12lb = F12*NewtonToLbConversion;

%% Part C
% Find the linear momentum of the system and the velocity 
% of the center of mass
% using L == linear momentum

L12 = m1*r1_dot + m2*r2_dot;
Vcm = L12 / (m1 + m2);

%% Part D
% What is the vector to the center of mass 10 seconds after these
% conditions?

Rcm_10s = Rcm + Vcm*10;

%% Part E
% If at another instant R2_dot is given as [0; 1.00016; 0], what is R1_dot?

% In this part, since momentum is a constant and conserved, I think I can
% solve for the new R1_dot
% m1R1_dot + m2R2_dot = p12

R2_dot_new = [0; 1.00016; 0];
R1_dot_new = (L12 - m2*R2_dot_new)/(m1);

%% Part F
% Calculate r12, r12_dot and mu
r12 = r2 - r1;
r12_dot = r2_dot - r1_dot;
mu = G*(m1+m2);

%% Part G
% Calculate a,e and Period of m2 relative to m1
% First I will find the specific energy for 
% r12_dot and r_12 distance, using mu

E = dot(r12_dot,r12_dot)/2 - mu/norm(r12); 

% Then using the specific energy = -mu/2a
a = -mu /(2*E) ;

h_norm = norm(cross(r12,r12_dot));
e = sqrt(1 - (h_norm)^2/(mu * a));
Period = 2*pi*sqrt(a^3/mu);

%% Part H
% How much would the smaller body's velocity have
% to increase at its current distance from the 
% larger body to reach escape velocity?

% At parabolic trajectory (needed to escape), the 
% energy of the system is 0. So V_esc = sqrt(2mu/r12)
% our delta-v will be Vesc - r12_dot

Vesc = sqrt(2*mu/norm(r12));
delta_v = Vesc - norm(r12_dot);



