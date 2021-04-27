clear, clc, close all;

% Driver Script for Module 11 Homework
% Dillon Allen

%% Problem 1

clear, clc, close all;

% Problem 1
% Find:
% 1. Initial State Error
% 2. Estimated State from the Least Squares Method
% 3. Final State Error
% 4. Discussion on improvement

load('HW11_Prob_1_data.mat');

p1initialStateError = norm(xVecTrue - xInitial);
H = [1 0; 1 1; 1 2; 1 4; 1 5];
x_hat = inv(H'*H)*H'*zVec;
p1finalStateError = norm(xVecTrue - x_hat);

% The values did not improve in this case

%% Problem 2
clear, clc, close all;

% Rolling Ball Weighted Least Squares
% Find:
% 1. Initial State Error
% 2. Estimated State from Least Squares (new data so not the same)
% 3. State Error from Least Squares
% 4. Estimated State from WLS
% 5. State error
% 6. Which was better?

load('HW11_Prob_2_data.mat');

initialStateError = norm(xVecTrue - xInitial);
H = [1 0; 1 1; 1 2; 1 3; 1 4];
x_hat_LS = inv(H'*H)*H'*zVec;
errorLS = norm(xVecTrue - x_hat_LS);

w = diag(sigmaVec.^2);

x_hat_WLS = inv(H'*inv(w)*H)*(H'*inv(w)*zVec);
errorWLS = norm(xVecTrue - x_hat_WLS);

% WLS was a better algorithm to use because it takes into account the
% weights per measurement

%% Problem 3

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

% Yes, it coverged pretty quickly

%% Problem 4
clear, clc, close all;

% Problem 4
% NLWLS AzEl Problem
% Use the provided functions for 2bodyECI, AzElECI and
% frame_ECEF2TangentENV
% .dat file did not have x_initial, so assuming x0 = x_initial
% P matrix is out of whack so I know there is an error here
% Final state error is of 10^7 power so also an issue

load('HW11_Prob_4_data.mat');

wVec = [1; 1;].* sigma_angle_rad;
W = eye(6) .* sigma_angle_rad^2;
Rinv = inv(W);

initStateError = norm(xVecTrue - x0);

x_star = x0;
opts1 = opts.stationOpts{1,1};
opts2 = opts.stationOpts{1,2};
opts3 = opts.stationOpts{1,3};
u = [];
v = [];
params = 0;

for i = 1 : length(tVec) - 1
    tspan = [tVec(i) tVec(i+1)];
    [t, x_star, Phi, Gamma_u, Gamma_v] = Provided_dynFxn_2Body_ECI(x0,u,v,tspan,params);
    [h1, H1, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, wVec, tspan(2), opts1);
    [h2, H2, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, wVec, tspan(2), opts2);
    [h3, H3, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, wVec, tspan(2), opts3);
    h = [h1; h2; h3];
    H = [H1; H2; H3];
    H = H*Phi;
    x_hat = x_star - inv(H' * Rinv * H)*(H'*Rinv*(zVecHist(i) - h));
    if(norm(x_hat - x_star) < 1e-6)
        disp('converge');
        break;
    else
        x_star = x_hat;
    end
end

P = inv(H'*Rinv*H);

finalError = norm(xVecTrue - x_star);

%% Function for truncating


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