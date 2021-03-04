clear, clc, close all;  

% Use r2 and v2 and calculate the six classical orbital elements. Confirm
% that these elements, except for theta, match the Starman example in
% Homework 3.

% Values from p1
r2 = [-135862465.167508, 62272677.7235982, -882257.467812311]'; % km
v2 = [-16.0245963820145 -29.1408936759632 -0.606602927396935]'; % km/s
mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

% Unit Vectors
I = [1 0 0]';
K = [ 0 0 1]';

v = norm(v2);
r = norm(r2);

E = (v^2)/2 - mu_S/r;
a = -mu_S/(2*E); % km
aAU = a * (1/1.498E8) % AU

hvec = cross(r2, v2);
evec = cross(v2,hvec)/mu_S - r2/r;
e = norm(evec);

i = acos(dot(hvec,K)/norm(hvec));

Nhat = cross(K,hvec)/norm(cross(K,hvec));
Omega = acos(dot(Nhat,I));
if Nhat(2) < 0
    Omega = 2*pi - Omega;
end

omega = acos(dot(Nhat,evec)/e);
if evec(3) < 0
   omega = 2*pi - omega; 
end

theta = acos(dot(evec,r2)/(e*r));
if dot(r2,v2) < 0
    theta = 2*pi - theta;
end

