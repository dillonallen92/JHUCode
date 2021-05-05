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



%% Plot
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