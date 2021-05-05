clear, clc, close all;

load('HW11_Prob_4corr.mat');

opts1 = opts.stationOpts{1,1};
opts2 = opts.stationOpts{1,2};
opts3 = opts.stationOpts{1,3};

initError = norm(xVecTrue - x0);

w = [0; 0];
R = eye(6) * sigma_angle_rad^2;
Rinv = inv(R); % W

% Following Batch processing flow chart, no prioris
lambda = 0; N = 0;
u = []; v = [];
param = [];
itNum = 0;
while 1
    % x_star = x0;
    for i = 1 : length(tVec)-1
        tspan = [tVec(1) tVec(i+1)];
        [t, x_star, Phi, Gamma_u, Gamma_v] = Provided_dynFxn_2Body_ECI(x0, u, v, tspan, param);
        [h1, H1, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, w, tspan(2), opts1);
        [h2, H2, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, w, tspan(2), opts2);
        [h3, H3, ~] = Provided_measFxnAzEl_ECI_HW11(x_star, w, tspan(2), opts3);
        h = [h1; h2; h3];
        Htilde = [H1; H2; H3];
        y(:,i) = zVecHist(:,i) - h;
        H = Htilde * Phi;
        lambda = lambda + transpose(H)*Rinv*H;
        N = N + transpose(H)*Rinv*y(:,i);
    end

    % Solve the normal equation
    x_hat = lambda\N
   
    if(norm(x_hat) < 1e-6)
        disp('converged');
        break;
    elseif itNum > 100
        disp('run time out');
        break;
    else
        x0 = x0 + x_hat;
        itNum = itNum + 1;
    end
end



