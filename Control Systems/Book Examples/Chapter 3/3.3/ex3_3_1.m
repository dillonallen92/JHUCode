% This code is a verification of the handwork done for example 3-3-1 in the
% book

% Once you have the transfer function foiled out on the numerator and
% denominator, we can use the following definitions

% numerator (5s + 3)
b = [5 3]; %numerator coefficients

% denominator (s^3 + 6s^2 + 11s + 6)
a = [1 6 11 6]; %denominator coefficients

% We can calculate the partial fraction expansion by using residue(b,a)
[r,p,k] = residue(b,a);

% The resulting values we get will be the following:
%   r -> the coefficients of the partial fraction
%   p -> the poles for the partial fractions
%   k -> polynomial (mostly 0 or constant [MATLAB DOC])

% Lets say we have the r,p,k values and we want to generate the polynomial
% coefficients, we can do the process in reverse

[b,a] = residue(r,p,k);
