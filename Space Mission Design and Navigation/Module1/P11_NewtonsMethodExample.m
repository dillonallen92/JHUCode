clear, clc, close all;
format long g
% Newton's Method Example from Notes
% Constant Variables
M = 0.25;
e = 0.1;

f = @(x) M - x + e*sin(x);
g = @(x) -1 + e*cos(x);

E = zeros(1001);
E(1) = M;

for i = 1:1000
   E(i+1) = E(i) - f(E(i))/g(E(i));
   if(f(E(i+1)) < 1E-15)
      fprintf('Eccentric Anomaly: %d\n', E(i+1));
      break;
   end
end