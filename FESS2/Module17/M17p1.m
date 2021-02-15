clear, clc, close all;

% Problem 1 part 1 of Module 17 Problem Set
% Given the Sun and Nadir Directions in inertial coords
Si = [0;0;1];
Ei = [1;0;0];

% Given Measurements of the same directions but in body frame
Sb = [0;1;0];
Eb = [0;0;1];

% Use the sun vector for the first body vector

TRIAD_Algo(Sb, Eb, Si, Ei)

%% Problem 2
Eb2 = [1/sqrt(2); 0; 1/sqrt(2)];

[attitude, axis, theta] = TRIAD_Algo(Sb, Eb2, Si, Ei)


%%