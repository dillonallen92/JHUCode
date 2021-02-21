clear, clc, close all;

% This script is to test function 2 to recreate the starman position and
% velocity vectors

a = 1.32489106722386;
e = 0.255915584525353;
omega = 3.09848024449946;
i = 0.0188055315710807; 
Omega = 5.53414135342721;
theta = 3.88465273236322;

%r is AU, v is AU / day
[r, v] = problem2Function(a,e,i,Omega,omega,theta)