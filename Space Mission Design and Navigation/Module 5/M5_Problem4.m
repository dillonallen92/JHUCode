clear, clc, close all;

% Solve the transfer with Lambert's Problem

% Info from problem 3
r1 = [-131386230.977293, 69971484.9501445, -718889.822774674]'; % km
r3 = [-139952788.024352, 54396927.5494364, -1043119.92803938]'; % km

% Sun gravitational parameter
mu_S = 1.327124400189E20 * (1/1000)^3; % km^3/s^2

% Ask the professor but I believe I can use my a from problem 2

r1norm = norm(r1); % km
r3norm = norm(r3); % km
a = 1.982008825664394e+08; % km (from problem 2)
c = norm(r3-r1); % km
s = (r1norm + r3norm + c)/2; % km
alpha = 2*asin(sqrt(s/(2*a))); % radians 
beta = 2*asin(sqrt((s-c)/(2*a))); % radians

%% Part A
% Find the minimum semi-major axis and the corresponding transfer time to
% go from r1 to r3

a_min = s/2; % km 
beta_min = 2*asin(sqrt((s-c)/(s))); % radians

tmin = sqrt(s^3 / (8*mu_S))*(pi - beta_min + sin(beta_min)); % s
minDays = tmin*(1/(3600*24)); % days

%% Part B
% Find v1 by solving Lambert's problem for a transfer from r1 to r3 using
% the transfer time found in Problem 3. Use Matlab's fzero function.

tProb3 = 5.347366092834404e+05; % s

theta = acos(dot(r1,r3)/(r1norm*r3norm)); % radians

[sol, fval, exitflag, output] = fzero(@(x) lamFun(x,tProb3,s,c,theta,tmin,mu_S), [a_min,10*a_min]);

alpha_vel = 2*asin(sqrt(s/(2*sol)));
beta_vel = 2*asin(sqrt((s-c)/(2*sol)));
A = sqrt(mu_S/(4*sol))*cot(alpha_vel/2); % km/s
B = sqrt(mu_S/(4*sol))*cot(beta_vel/2); % km/s 
uc = (r3-r1)/norm(r3-r1); 
u1 = r1/r1norm;

v1 = (B+A)*uc + (B-A)*u1; % km/s

%% Part C
% Show that the resulting value of semi-major axis matches that of Starman
aSM = 198200882.566439; % km (StarMan)

if aSM - sol < 1E-5
    disp("close enough")
end


function val = lamFun(x,tof,s,c,theta,tmin, mu_S)

    alpha_0 = 2*asin(sqrt(s/(2*x)));
    beta_0 = 2*asin(sqrt((s-c)/(2*x)));

    if (0 <= theta) && (theta <= pi)
        beta = beta_0;
    else
        beta = - beta_0;
    end
    
    if tof <= tmin
        alpha = alpha_0;
    else
        alpha = 2*pi - alpha_0;
    end
    
    val = sqrt(mu_S/x^3)*tof -  (alpha - beta - (sin(alpha) - sin(beta)));
end

