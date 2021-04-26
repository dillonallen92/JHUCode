% Problem 3
% NLS Range

clear, clc, close all;

load('HW11_Prob_3_data.mat');

opts1.r_IVec = [r1_IVec(1); r1_IVec(2);];
opts2.r_IVec = [r2_IVec(1); r2_IVec(2);];
xVec = zVec;
t = 0;
w = sigmaMeas;
W = eye(2) * sigmaMeas^2;
Rinv = inv(W);

% Truncate xVecTrue to remove z,vx,vy,vz dependance
xVecTrue_trunc = [xVecTrue(1); xVecTrue(2)];

% Initiate x_star
x_star = x0;
j = 1;

% Initial State Error
initError = norm(xVecTrue_trunc - x0);

% Loop til convergence
for  j  = 1 : 10
    [h1, H1, ~] = Provided_measFxnRange(x_star, w, t, opts1);
    [h2, H2, ~] = Provided_measFxnRange(x_star, w, t, opts2);
    h = [h1; h2];
    H = [H1; H2];
    x_hat = x_star - inv(H'*Rinv *H)*H'* Rinv* (h - zVec);
    xVecHist(:,j) = x_hat;
    errHist(j) = norm(x_hat - xVecTrue_trunc);
    x_star = x_hat
end

% P Matrix
P = inv(H'*Rinv*H);

% Final State Error
finalErr = errHist(length(errHist));

% Plot of error vs iteration count
itVec = 1:1:length(errHist);
plot(errHist, itVec);
xlabel("Iteration #");
ylabel("Estimated State Error");
title("Estimated Error vs Iteration Count");

% Function to test truncating
function [h,HMat,id] = Provided_measFxnRange(xVec,w,t,opts)
    % Parse inputs:
    rVec   = xVec(1:2,1);
    % n      = size(xVec,1);
    r_IVec = opts.r_IVec;
    % Compute rho:
    dr     = rVec - r_IVec;
    rho    = sqrt( dr'*dr );
    % Compute H:
    HMat   = [ (rVec(1)-r_IVec(1))/rho, (rVec(2)-r_IVec(2))/rho];

    % If truth-model, then add noise:
    h = rho;
    if ~isempty(w)
        h = h + w;
    end

    % ID: (not needed yet)
    id = [];

    return
end