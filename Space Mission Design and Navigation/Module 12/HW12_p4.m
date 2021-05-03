clear, clc, close all;

load('HW12_Prob_4.mat');

% Ideal Range Extended Kalman Filter Problem
% Use the provided dynamics and measurement models provided
% dynFxn_2Body_ECI
% measFxnRange_ECI

initError = norm(x0_True - xStar_0Vec);

xVecHist = zeros(length(xStar_0Vec), length(tVec));
xVecHist(:,1) = xStar_0Vec;
x_hat(:,1) = xStar_0Vec;
u = []; v = []; params = [];
P{1,1} = PBarMat0;
w = 0;

for i = 1 : length(tVec) - 1
    
   % Set the timespan
   tspan = [tVec(i), tVec(i+1)];
   
   % Compute the dynamics
   [t, xVecHist(:,i+1), Phi, Gamma_u, Gamma_v] = Provided_dynFxn_2Body_ECI(x_hat(:,i), ...
                                                                           u,v,tspan,params);
   
   % Prop the covariance and handle linearization points
   Pbar = Phi * P{1,i} * transpose(Phi) + Gamma_u * QMat * transpose(Gamma_u);
   [h, H, ~] = Provided_measFxnRange_ECI(xVecHist(:,i+1),w,tspan(2),measOpts);
   y(i) = zVecHist(i) - h;
   
   % Execute the Kalman Filter
   K = Pbar*transpose(H)*inv(H*Pbar*transpose(H) + RMat);
   x_hat(:,i+1) = xVecHist(:,i+1) + K*y(i);
   kh = K*H;
   P{1,i+1} = (eye(size(kh)) - kh)*Pbar;
end

for i = 1 : length(xTrueVecHist)
   errorVec(:,i) = xTrueVecHist(:,i) - x_hat(:,i); 
   errorRes(i) = norm(errorVec(:,i));
   xVecPosStd(i) = sqrt(P{1,i}(1,1));
   xVecNegStd(i) = -sqrt(P{1,i}(1,1));
end

plot(tVec, errorVec(1,:)); hold on;
plot(tVec, xVecPosStd);
plot(tVec, xVecNegStd);
xlabel("time (s)");
ylabel("Xpos (m)");
hold off;

