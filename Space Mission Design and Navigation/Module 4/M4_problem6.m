clear, clc, close all;

% Problem 1 in Tapley book

% Initial Values
x0 = 1.5;
y0 = 10.0;
xdot0 = 2.2;
ydot0 = 0.5;
g = 0.3;
xs = 1.0;
ys = 1.0;

% range observations
p0 = 7.0;
p1 = 8.00390597;
p2 = 8.94427191;
p3 = 9.801147892;
p4 = 10.630145813;

c = [p0 p1 p2 p3 p4]';

xvec0 = [x0 y0 xdot0 ydot0 g]';
tol = 1E-5;

xp = NewtonScheme(c, xvec0, xs, ys, tol);

disp("-- Actual Initial Conditions --")
fprintf("x0: %d\n", xp(1));
fprintf("y0: %d\n", xp(2));
fprintf("xdot0: %d\n", xp(3));
fprintf("ydot0: %d\n", xp(4));
fprintf("g0: %d\n", xp(5));

% Functions

% X component helper function
function val = xcomp(xvec, t, xs)
    val = xvec(1) - xs + xvec(3)*t;
end

% Y component helper function
function val = ycomp(xvec,t,ys)
    val = xvec(2) - ys + xvec(4)*t - .5*xvec(5)*t^2;
end

% Range Equation
function val = rho(xvec, t, xs, ys)
    val = sqrt(xcomp(xvec,t, xs)^2 + ycomp(xvec, t,ys)^2);
end

% Initial G vector function
function val = G(xvec,xs,ys)
    val = zeros(5,1);
    for i = 1 : 5
        t = i - 1;
        val(i) = rho(xvec, t, xs, ys);
    end
end

% Matrix of partial derivatives
function val = pJ(xvec, xs, ys)
    val = zeros(5); % 5 x 5 matrix
    for i = 1 : 5
       t = i - 1;
       val(i,1) = xcomp(xvec, t, xs) / rho(xvec, t, xs, ys) ;
       val(i,2) = ycomp(xvec, t, ys) / rho(xvec, t, xs, ys) ;
       val(i,3) = (t * xcomp(xvec, t, xs)) / rho(xvec, t, xs, ys);
       val(i,4) = (t * ycomp(xvec, t, ys)) / rho(xvec, t, xs, ys);
       val(i,5) = (-t^2*ycomp(xvec, t, ys))/(2*rho(xvec, t, xs, ys));
    end
end

% Newtons Method Function
function xp = NewtonScheme(c, xvec, xs, ys, tol)
    while 1
       J = c - G(xvec, xs, ys);
       xp = xvec + pJ(xvec, xs, ys)\J; % My pJ is - pG/pX0 so there is a +
       if norm(xp - xvec) < tol
           break;
       end
       xvec = xp;
    end
end