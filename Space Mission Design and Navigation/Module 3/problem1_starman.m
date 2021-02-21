clear, clc, close all;
format long g
% This script is to test the Problem1Function algorithm that 
% is designed to calculate the six orbital elements for the 
% starman example

r = [1.523679990458238; -7.50002657546034e-2; 1.848063114121444e-2]; % AU
v = [-2.058658387423061e-3; 1.266141033930193e-2; 1.480290066794905e-4]; %AU/Day

[a_au,e_norm,i,omega,Omega,theta] = problem1Function(r,v);