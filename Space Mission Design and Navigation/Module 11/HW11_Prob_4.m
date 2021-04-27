clear, clc, close all;

% Problem 4
% NLWLS AzEl Problem
% Use the provided functions for 2bodyECI, AzElECI and
% frame_ECEF2TangentENV
% .dat file did not have x_initial, so assuming x0 = x_initial

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