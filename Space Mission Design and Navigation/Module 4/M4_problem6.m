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

while 1
   J = c - G(xvec0, xs, ys);
   xp = xvec0 + pJ(xvec0, xs, ys)\J; % My pJ is - pG/pX0 so there is a +
   if norm(xp - xvec0) < tol
       disp('Real initial values:');
       xp
       break;
   end
   xvec0 = xp;
end

function val = xcomp(xvec, t, xs)
    val = xvec(1) - xs + xvec(3)*t;
end

function val = ycomp(xvec,t,ys)
    val = xvec(2) - ys + xvec(4)*t - .5*xvec(5)*t^2;
end

function val = rho(xvec, t, xs, ys)
    val = sqrt(xcomp(xvec,t, xs)^2 + ycomp(xvec, t,ys)^2);
end

function val = G(xvec,xs,ys)
    val = zeros(5,1);
    for i = 1 : 5
        t = i - 1;
        val(i) = rho(xvec, t, xs, ys);
    end
end

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