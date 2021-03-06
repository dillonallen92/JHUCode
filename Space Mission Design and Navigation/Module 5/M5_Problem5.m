clear, clc, close all;

% Constants
mu_S = 1.327124400189E20 * (1/1000)^3; 
tspan=linspace(0,5.347366092834404e+05,100); % s

r1 = [-131386230.977293, 69971484.9501445, -718889.822774674]'; % km
v1 = [-17.453946908042990;-28.433493882607780;-0.615165025022269]; % km /s

x0vec = [r1; v1];

integOptions=odeset('abstol',1e-10,'reltol',1e-10);
[t,x]=ode113(@(t,x)twoBodyCart(t,x,mu_S), tspan, x0vec, integOptions);

% Problem 4 Prop
function dx = twoBodyCart(~, x, mu_in)
    dx = zeros(6,1);
    dx(1) = x(4);
    dx(2) = x(5);
    dx(3) = x(6);
    dx(4) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(1);
    dx(5) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(2);
    dx(6) = (-mu_in / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2))*x(3);
end
