clear, clc, close all;

%% Load the data
load('HW8_Prob3.mat');

%% Part 1
% Plot the electron density vs altitude, with altitude in linear units
% on the y-axis and log(density) units on the x-axis.
z = @(h) (h - (hmax_m/1000))/((hscale_m/1000)); % converted to km
rho = @(z) pmax_epm3 * exp(1 - z - exp(-z)); % e- per m^3
alt_max = 2000;
height = zeros(1,alt_max);
density_mat = zeros(1,alt_max);

for i = 0 : alt_max
    height(i+1) = i;
    density_mat(i+1) = rho(z(i));
end

plot(log10(density_mat), height);
xlabel("log(\rho)");
ylabel("Altitude (km)");
title("e^{-} Density Profile");
%% Part 2
% Integrate the TEC vertically through the profile
VTEC = trapz(height.*1000, density_mat);

%% Part 3
% What is the resulting signal delay at the GPS L1 frequency
% f = 1.57542 GHz due to this electron content?
% Note: use delT = 40.3*TEC/(c*f^2)

c = 3e8; % m /s
f = 1.57542e9; % 1/s

delT = 40.3*VTEC/(c*f^2);