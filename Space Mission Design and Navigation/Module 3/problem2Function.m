function [r,v] = problem2Function(a,e,i,Omega,omega,theta)

    % This function takes all six classical orbital elements
    % and determines the position and velocity vectors

    % Constants
    mu = 132712440041.94; % km^3/s^2
    KmInAU = 149597870.7; % AU/ km
    secondsInDay = 24*3600; % days / second
    
    % Using Prussing as a guide
    a = a * KmInAU;
    r_norm = (a*(1-e^2))/(1+e*cos(theta));
    theta = theta + omega; % theta in prussing is theta + omega for us
    
    r_x = r_norm*(cos(Omega)*cos(theta) - sin(Omega)*sin(theta)*cos(i));
    r_y = r_norm*(sin(Omega)*cos(theta) + cos(Omega)*sin(theta)*cos(i));
    r_z = r_norm*sin(theta)*sin(i);
    
    r = [r_x; r_y; r_z];
    r = r / KmInAU;
    
    % Velocity
    h = sqrt(mu*a*(1-e^2));
    sinTerms = sin(theta) + e*sin(omega);
    cosTerms = cos(theta) + e*cos(omega);
    muh = mu / h;
    
    v_x = -muh * (cos(Omega) * sinTerms + sin(Omega)*cosTerms*cos(i));
    v_y = muh * (sin(Omega)*sinTerms - cos(Omega)*cosTerms*cos(i));
    v_z = muh * cosTerms * sin(i);
    
    v = [v_x; v_y; v_z];
    
    v = v * secondsInDay / KmInAU;
    
    
end