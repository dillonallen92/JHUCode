function [r,v] = problem2Function(a,e,i,Omega,omega,theta)

    % This function takes all six classical orbital elements
    % and determines the position and velocity vectors

    % Constants
    mu = 132712440041.94; % km^3/s^2
    auToKm = 149597870.7; % km
    dayToSeconds = 1/(24*3600);

    % Using Prussing as a guide

    theta_epoch = omega + theta;
    r_norm = (a*(1-e^2))/(1+e*cos(theta_epoch));

    
end