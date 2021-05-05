% HW12 Driver Script for all problems
% Dillon Allen
% 5 / 4 / 2021

%% Problem 1
% -- Rolling Ball problem --

clear, clc, close all;

load('HW12_prob_1.mat');

xStarVecHist = zeros(size(xTrueVecHist));
x_hat = zeros(size(xTrueVecHist));
xStarVecHist(:,1) = xStar_0Vec;
x_hat(:,1) = xStar_0Vec;

Pmat{1,1} = PBarMat0;
xPos3Std(1) = 3 * sqrt(Pmat{1,1}(1,1));
xNeg3Std(1) = -3 * sqrt(Pmat{1,1}(1,1));

% Initial State Error
errHist(:,1) = xStar_0Vec - xTrue0;

for i = 1 : length(tVec) - 1
   % Model formulas
   dt = tVec(i);
   H = [1 0];
   Phi = [1 dt; 0 1];
   Gamma = [.5 * dt^2; dt];
   
   % Propagate the state
   xStarVecHist(:,i+1) = Phi * x_hat(:,1);
   
   % Propagate the covariance
   Pbar = Phi * Pmat{1,i} * transpose(Phi) + Gamma * QMat * transpose(Gamma);
   y(i) = zVecHist(i) - H * xStarVecHist(:,i+1);
   
   % Execute the Kalman Filter
   K = Pbar * transpose(H) * inv(H * Pbar * transpose(H) + RMat);
   x_hat(:,i+1) = xStarVecHist(:,i+1) + K * y(i);
   kh = K*H;
   Pmat{1,i+1} = (eye(size(kh)) - kh) * Pbar;
   
   errHist(:,i+1) =  x_hat(:,i+1) - xTrueVecHist(:,i+1);
   xPos3Std(i+1) = 3 * sqrt(Pmat{1,i+1}(1,1));
   xNeg3Std(i+1) = -3 * sqrt(Pmat{1,i+1}(1,1));
   vPos3Std(i+1) = 3 * sqrt(Pmat{1,i+1}(2,2));
   vNeg3Std(i+1) = -3 * sqrt(Pmat{1,i+1}(2,2));
end



% Plot
for i = 1 : 2
    
   subplot(1,2,i)
   plot(tVec, errHist(i,:)); hold on;
   if i == 1
       plot(tVec, xPos3Std);
       plot(tVec, xNeg3Std);
       xlabel("time (s)");
       ylabel("Position Error (m)");
   else
       plot(tVec, vPos3Std);
       plot(tVec, vNeg3Std);
       xlabel("time (s)");
       ylabel("Velocity Error (m/s)");
   end
end

% Compare to propagating the state
% Looking at x_hat, xTrueVecHist and xStarVecHist, you can see that x_hat
% and xStarVecHist approach the same final values while xStarVecHist
% overshoots. The estimate is therefore an improvement from just doing
% propagation.

%% Problem 2
% -- Theory Question --
% The QMat your boss gives you will be higher because the disturbances at
% the lower altitude should be higher. Therefore, this is not the correct
% matrix for your model and the values you need should be smaller.

%% Problem 3
% -- TDOA/FDOA Model -- 

clear, clc, close all;

load('HW12_Prob_3.mat');

xVecHist = zeros(length(xStar_0Vec), length(tVec));
xVecHist(:,1) = xStar_0Vec;
x_hat(:,1) = xStar_0Vec;
u = []; v = []; params = [];
P{1,1} = PBarMat0;
w = 0;

% Initial State Error
initErr = xStar_0Vec - x0_True;

for i = 1 : length(tVec) - 1
   % create the tspan vec
   tspan = [tVec(i), tVec(i+1)];
   
   % Propagate the state
   [t, xVecHist(:,i+1), Phi, Gamma_u, Gamma_v] = Provided_dynFxn_2Body_ECI(x_hat(:,i),u,v,tspan,params);
   
   % Compute the measurement model
   [h, H, ~] = Provided_measFxn_TDOA_and_FDOA_ECEFStat(xVecHist(:,i+1),w,tspan(2), measOpts);
   y(:,i) = zVecHist(:,i) - h;
   
   % Propagate Covariance
   Pbar = Phi * P{1,i} * transpose(Phi) + Gamma_u * QMat * transpose(Gamma_u);
   
   % Execute the Kalman Filter Equations
   K = Pbar * transpose(H) * inv(H * Pbar * transpose(H) + RMat);
   x_hat(:,i+1) = xVecHist(:,i+1) + K*y(:,i);
   kh = K*H;
   P{1,i+1} = (eye(size(kh)) - kh)* Pbar;
end

% Plotting stuff

for i = 1 : length(xTrueVecHist)
   errorVec(:,i) = xTrueVecHist(:,i) - x_hat(:,i); 
   xVecPosStd(i) = sqrt(P{1,i}(1,1));
   xVecNegStd(i) = -sqrt(P{1,i}(1,1));
   yVecPosStd(i) = sqrt(P{1,i}(2,2));
   yVecNegStd(i) = -sqrt(P{1,i}(2,2));
   zVecPosStd(i) = sqrt(P{1,i}(3,3));
   zVecNegStd(i) = -sqrt(P{1,i}(3,3));
   xdotVecPosStd(i) = sqrt(P{1,i}(4,4));
   xdotVecNegStd(i) = -sqrt(P{1,i}(4,4));
   ydotVecPosStd(i) = sqrt(P{1,i}(5,5));
   ydotVecNegStd(i) = -sqrt(P{1,i}(5,5));
   zdotVecPosStd(i) = sqrt(P{1,i}(6,6));
   zdotVecNegStd(i) = -sqrt(P{1,i}(6,6));
end


subplot(2,3,1)
plot(tVec, errorVec(1,:)); hold on;
plot(tVec, xVecPosStd);
plot(tVec, xVecNegStd);
xlabel("time (s)");
ylabel("Xpos (m)");
legend('Xpos', '+/- 1-\sigma');
hold off;

subplot(2,3,2)
plot(tVec, errorVec(2,:)); hold on;
plot(tVec, yVecPosStd);
plot(tVec, yVecNegStd);
xlabel("time (s)");
ylabel("Ypos (m)");
legend('Ypos', '+/- 1-\sigma');

subplot(2,3,3)
plot(tVec, errorVec(3,:)); hold on;
plot(tVec, zVecPosStd);
plot(tVec, zVecNegStd);
xlabel("time (s)");
ylabel("Zpos (m)");
legend('Zpos', '+/- 1-\sigma');

subplot(2,3,4)
plot(tVec, errorVec(4,:)); hold on;
plot(tVec, xdotVecPosStd);
plot(tVec, xdotVecNegStd);
xlabel("time (s)");
ylabel("xVel (m/s)");
legend('xVel', '+/- 1-\sigma');

subplot(2,3,5)
plot(tVec, errorVec(5,:)); hold on;
plot(tVec, ydotVecPosStd);
plot(tVec, ydotVecNegStd);
xlabel("time (s)");
ylabel("yVel (m/s)");
legend('yVel', '+/- 1-\sigma');

subplot(2,3,6)
plot(tVec, errorVec(6,:)); hold on;
plot(tVec, zdotVecPosStd);
plot(tVec, zdotVecNegStd);
xlabel("time (s)");
ylabel("zVel (m/s)");
legend('zVel', '+/- 1-\sigma');

sgtitle('Position and Velocity errors over time');

%% Problem 4
% -- Ideal Range -- 
% Ideal Range Extended Kalman Filter Problem
% Use the provided dynamics and measurement models provided
% dynFxn_2Body_ECI
% measFxnRange_ECI

clear, clc, close all;

load('HW12_Prob_4.mat');

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
   
   % Chi-Square Stuff
   
   % Innovation Eq
   q_meas(i) = transpose(y(i)) * inv(H * Pbar * transpose(H) + RMat) * y((i));
   
   % State Error
   q_error(i) = transpose(x_hat(:,i) - xTrueVecHist(:,i)) * inv(P{1,i}) * (x_hat(:,i) - xTrueVecHist(:,i));
end

% Plotting 

for i = 1 : length(xTrueVecHist)
   errorVec(:,i) = xTrueVecHist(:,i) - x_hat(:,i); 
   xVecPosStd(i) = sqrt(P{1,i}(1,1));
   xVecNegStd(i) = -sqrt(P{1,i}(1,1));
   yVecPosStd(i) = sqrt(P{1,i}(2,2));
   yVecNegStd(i) = -sqrt(P{1,i}(2,2));
   zVecPosStd(i) = sqrt(P{1,i}(3,3));
   zVecNegStd(i) = -sqrt(P{1,i}(3,3));
   xdotVecPosStd(i) = sqrt(P{1,i}(4,4));
   xdotVecNegStd(i) = -sqrt(P{1,i}(4,4));
   ydotVecPosStd(i) = sqrt(P{1,i}(5,5));
   ydotVecNegStd(i) = -sqrt(P{1,i}(5,5));
   zdotVecPosStd(i) = sqrt(P{1,i}(6,6));
   zdotVecNegStd(i) = -sqrt(P{1,i}(6,6));
end


subplot(2,3,1)
plot(tVec, errorVec(1,:)); hold on;
plot(tVec, xVecPosStd);
plot(tVec, xVecNegStd);
xlabel("time (s)");
ylabel("Xpos (m)");
legend('Xpos', '+/- 1 \sigma');
hold off;

subplot(2,3,2)
plot(tVec, errorVec(2,:)); hold on;
plot(tVec, yVecPosStd);
plot(tVec, yVecNegStd);
xlabel("time (s)");
ylabel("Ypos (m)");
legend('Ypos', '+/- 1 \sigma');

subplot(2,3,3)
plot(tVec, errorVec(3,:)); hold on;
plot(tVec, zVecPosStd);
plot(tVec, zVecNegStd);
xlabel("time (s)");
ylabel("Zpos (m)");
legend('Zpos', '+/- 1 \sigma');

subplot(2,3,4)
plot(tVec, errorVec(4,:)); hold on;
plot(tVec, xdotVecPosStd);
plot(tVec, xdotVecNegStd);
xlabel("time (s)");
ylabel("xVel (m/s)");
legend('xVel', '+/- 1 \sigma');

subplot(2,3,5)
plot(tVec, errorVec(5,:)); hold on;
plot(tVec, ydotVecPosStd);
plot(tVec, ydotVecNegStd);
xlabel("time (s)");
ylabel("yVel (m/s)");
legend('yVel', '+/- 1 \sigma');

subplot(2,3,6)
plot(tVec, errorVec(6,:)); hold on;
plot(tVec, zdotVecPosStd);
plot(tVec, zdotVecNegStd);
xlabel("time (s)");
ylabel("zVel (m/s)");
legend('zVel', '+/- 1 \sigma');

sgtitle('Position and Velocity errors over time');

% Chi-Square stuff (?)

y_meas = chi2pdf(q_meas, 3);
y_error = chi2pdf(q_error, 3);

subplot(1,2,1);
scatter(q_meas, y_meas);
xlabel("q_meas");
ylabel("q_meas pdf");

subplot(1,2,2);
scatter(q_error, y_error);
xlabel("q_error");
ylabel("q_error pdf");

% This looks properly distributed but honestly I am not exactly sure what I
% am supposed to be looking for with this besides seeing the curve I
% wanted.
