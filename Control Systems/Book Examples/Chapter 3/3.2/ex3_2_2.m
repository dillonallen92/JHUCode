% Use the "zeros poles gain" (zpk) function to create the transfer function
% polynomial model

G = zpk([-2], [0 -1 -3 -3], 10)

% Now use tf() to change G into a transfer function in Matlab
Gp = tf(G)

% You can find the zeros by using the command zero() and poles by using the
% command pole()

zeros = zero(Gp)
poles = pole(Gp)