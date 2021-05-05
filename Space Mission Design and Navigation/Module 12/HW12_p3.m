clear, clc, close all;

load('HW12_Prob_3.mat');

xVecHist = zeros(length(xStar_0Vec), length(tVec));
xVecHist(:,1) = xStar_0Vec;
x_hat(:,1) = xStar_0Vec;
u = []; v = []; params = [];
P{1,1} = PBarMat0;
w = 0;

% Initial State Error
initErr = norm(xStar_0Vec - x0_True);

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


