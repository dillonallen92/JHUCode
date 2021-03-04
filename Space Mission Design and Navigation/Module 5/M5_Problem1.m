clear, clc, close all;

% Position vectors about the sun
r1 = [-131386230.977293, 69971484.9501445, -718889.822774674]'; % km
r2 = [-135862465.167508, 62272677.7235982, -882257.467812311]'; % km
r3 = [-139952788.024352, 54396927.5494364, -1043119.92803938]'; % km

% Sun gravitational parameter
mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

%% Part A
% Use Gibb's method to determine p, e, h and v2

r1norm = norm(r1); % km
r2norm = norm(r2); % km
r3norm = norm(r3); % km

Nvec = r1norm*cross(r2,r3) + r2norm*cross(r3,r1) + r3norm*cross(r1,r2);
Dvec = cross(r1,r2) + cross(r2, r3) + cross(r3,r1);

N = norm(Nvec);
D = norm(Dvec);

p = N/D; % km

Svec = (r2norm - r3norm)*r1 + (r3norm - r1norm)*r2 + (r1norm - r2norm)*r3;
S = norm(Svec);

e = S/D;

h = sqrt(N*mu_S/D);

v2vec = (1/r2norm)*sqrt(mu_S/(N*D))*cross(Dvec,r2) + sqrt(mu_S/(N*D))*Svec; % km/s

%% Part B
% Additionally, use the equation for velocity to find v1 and v3

v1vec = (1/r1norm)*sqrt(mu_S/(N*D))*cross(Dvec,r1) + sqrt(mu_S/(N*D))*Svec; % km/s 
v3vec = (1/r3norm)*sqrt(mu_S/(N*D))*cross(Dvec,r3) + sqrt(mu_S/(N*D))*Svec; % km/s
