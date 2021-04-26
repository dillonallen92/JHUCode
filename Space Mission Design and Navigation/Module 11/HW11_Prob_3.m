% Problem 3
% NLWLS Range

clear, clc, close all;

load('HW11_Prob_3_data.mat');

opts1.r_IVec = [r1_IVec(1); r1_IVec(2);];
opts2.r_IVec = [r2_IVec(1); r2_IVec(2);];
xVec = zVec;
t = 0;
w = 0;

% Initiate x_star
x_star = x0;

% Loop til convergence
while 1
    [h1, H1, ~] = Provided_measFxnRange(x_star, w, t, opts1);
    [h2, H2, ~] = Provided_measFxnRange(x_star, w, t, opts2);
    h = [h1; h2];
    H = [H1; H2];
    x_hat = x_star - (H'*H)*H'*(h - zVec);
    if (norm(x_hat - x_star) < 1e-6)
        disp('converge');
        xVec_converge = x_hat;
        break;
    else
        x_star = x_hat
    end
end

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