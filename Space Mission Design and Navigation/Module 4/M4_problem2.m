clear, clc, close all

% Write matlab code that integrates the cartesian equations of motion

% constant
muE=3.986e5; % km^3/s^2
tspan=linspace(0,36*3600,100); % s

%% Part A
% Convert the initial conditions given for the circular orbit in the
% lecture into cartesian form and integrate over a timespan of 36 hours.
% Verify that your results match the result in the lecture.

% I am assuming z = 0 ? not sure what to do with that to be honest
% x = r*cos(theta);
% y = r*sin(theta);
% z = 0;
% xdot = r_dot * cos(theta) - r*sin(theta)*theta_dot;
% ydot = r_dot*sin(theta) + r*cos(theta) * theta_dot;
% zdot = 0;

% Initial Circle polar matrix from lecture
x0Circ=[42241.0800678832;0;0;7.27220521664305e-05];

% grab values from circle matrix
r = x0Circ(1);
theta = x0Circ(2);
r_dot = x0Circ(3);
theta_dot = x0Circ(4);

% polar to cartesian conversion
x = r*cos(theta); % km 
y = r*sin(theta); % km
z = 0; % km
xdot = r_dot * cos(theta) - r*sin(theta)*theta_dot; % km / s
ydot = r_dot*sin(theta) + r*cos(theta) * theta_dot; % km / s
zdot = 0; % km / s

x0CircCart = [x y z xdot ydot zdot]';

integOptions=odeset('abstol',1e-10,'reltol',1e-10);
[t,x]=ode113(@(t,x)twoBodyCart(t,x,muE),tspan,x0CircCart,integOptions);

plot(x(:,1), x(:,2));
xlabel(" x (km) " );
ylabel(" y (km) " );
title("Circular Plot");
axis equal


%% Part B
% Ellipse time

x0Ellipt=[21120.5400339416;0;0;0.000251916578365864]; % (km,rad,km/s,rad/s)

r = x0Ellipt(1);
theta = x0Ellipt(2);
r_dot = x0Ellipt(3);
theta_dot = x0Ellipt(4);

x = r*cos(theta); % km 
y = r*sin(theta); % km
z = 0; % km
xdot = r_dot * cos(theta) - r*sin(theta)*theta_dot; % km / s
ydot = r_dot*sin(theta) + r*cos(theta) * theta_dot; % km / s
zdot = 0; % km / s

x0ElliptCart = [x y z xdot ydot zdot]';

integOptions=odeset('abstol',1e-10,'reltol',1e-10);
[t,x]=ode113(@(t,x)twoBodyCart(t,x,muE),tspan,x0ElliptCart,integOptions);

plot(x(:,1), x(:,2));
xlabel(" x (km) " );
ylabel(" y (km) " );
title("Elliptical Plot");
axis equal

%% Function

function dx = twoBodyCart(~, x, mu_in)
    dx = zeros(6,1);
    dx(1) = x(4);
    dx(2) = x(5);
    dx(3) = x(6);
    dx(4) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(1);
    dx(5) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(2);
    dx(6) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(3);
end
