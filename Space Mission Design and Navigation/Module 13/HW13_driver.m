% Main script for all the work for HW13
% Problem 1 and 3 (reduced problem set)

clear, clc, close all;

%% Problem 1

load('HW13_Prob_1.mat');

% Rolling Ball Problem
% Kalman Filter and Smoother
% 1. Initial State Error
% 2. (EKF) Estimated States and Covariance Matrix time-histories
% 3. (RTS) Estimated States and Covariance Matrix time-histories 
% 4. Plot each state with +/- 3 sigma
% 5. Did the estimates improve?

% Variable setup
Pcell{1,1} = PBarMat0;
Pcell_smooth = cell(1,length(tVec)-1);
Pbarcell{1,1} = PBarMat0;
Phi_hist{1,1} = eye(2);
Smat = cell(1,length(tVec)-1);

x_hat(:,1) = xStar_0Vec;
x_hat_smooth = zeros(2,length(tVec)-1);
xVecHist(:,1) = xStar_0Vec;

% EKF
for i = 1 : length(tVec)-1
   tspan = tVec(1):tVec(i);
   dt = tVec(i);
   H = [1 0];
   Gamma = [.5*dt^2; dt];
   Phi = [ 1 dt; 0 1];
   Phi_hist{1,i+1} = Phi;
   xVecHist(:,i+1) = Phi * xVecHist(:,1);
   
   Pbar = Phi * Pcell{1,i} * transpose(Phi) + Gamma*QMat*transpose(Gamma);
   Pbarcell{1,i+1} = Pbar;
   yVec = zVecHist(i) - H*xVecHist(:,i+1);
   
   K = Pbar * transpose(H) * inv(H*Pbar*transpose(H) + RMat);
   x_hat(:,i+1) = xVecHist(:,i+1) + K*yVec;
   kh = K*H;
   Pcell{1,i+1} = (eye(size(kh)) - kh)*Pbar;
   
end

Pcell_smooth{1,length(tVec)} = Pcell{1,length(tVec)};
x_hat_smooth(:,length(tVec)) = x_hat(:,length(tVec));

% Smoother
for ell = length(tVec) : -1 : 2
    Smat{1,ell-1} = Pcell{1,ell-1} * transpose(Phi_hist{1,ell}) * inv(Pbarcell{1,ell});
    Pcell_smooth{1,ell-1} = Pcell{1,ell-1} + Smat{1,ell-1} * ...
                            (Pcell_smooth{1,ell} - Pcell{1,ell}) * transpose(Smat{1, ell-1});
    deltaX = x_hat_smooth(:,ell) - xVecHist(:,ell);
    x_hat_smooth(:,ell-1) = x_hat(:,ell-1) + Smat{1,ell-1} * deltaX;
end

% Error and Std data

for i = 1 : length(xTrueVecHist)
   % Error Matrices
   errHist(:,i) =  x_hat(:,i) - xTrueVecHist(:,i);
   errHistSmoothed(:,i) = x_hat_smooth(:,i) - xTrueVecHist(:,i);
   
   % Normal Covariances
   xPos3Std(i) = 3 * sqrt(Pcell{1,i}(1,1));
   xNeg3Std(i) = -3 * sqrt(Pcell{1,i}(1,1));
   vPos3Std(i) = 3 * sqrt(Pcell{1,i}(2,2));
   vNeg3Std(i) = -3 * sqrt(Pcell{1,i}(2,2)); 
   
   % Smoothed Covariances
   xPos3StdSmoothed(i) = 3 * sqrt(Pcell_smooth{1,i}(1,1));
   xNeg3StdSmoothed(i) = -3 * sqrt(Pcell_smooth{1,i}(1,1));
   vPos3StdSmoothed(i) = 3 * sqrt(Pcell_smooth{1,i}(2,2));
   vNeg3StdSmoothed(i) = -3 * sqrt(Pcell_smooth{1,i}(2,2));
end

% Plots

% X Vec
subplot(211);
plot(tVec, errHist(1,:), 'k-'); hold on; 
plot(tVec, errHistSmoothed(1,:), 'k--', 'LineWidth', 1)
plot(tVec, xPos3Std, 'b');
plot(tVec, xNeg3Std, 'b');
plot(tVec, xPos3StdSmoothed, 'r--','LineWidth', 0.8);
plot(tVec, xNeg3StdSmoothed, 'r--', 'LineWidth', 0.8);
xlabel("time (s)");
ylabel("X error");
lgnd = legend('unsmoothed', 'smoothed', ' \pm 3-\sigma unsmoothed', ...
       '\pm 3-\sigma unsmoothed', '\pm 3-\sigma smoothed', '\pm 3-\sigma smoothed');
lgnd.Location = 'northeast';
lgnd.BoxFace.ColorType='truecoloralpha'
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.55]')
lgnd.NumColumns = 2;
lgnd.Orientation = 'horizontal';
hold off;

% Xdot Vec
subplot(212)
plot(tVec, errHist(2,:), 'k-'); hold on; 
plot(tVec, errHistSmoothed(2,:), 'k--', 'LineWidth', 1)
plot(tVec, vPos3Std, 'b');
plot(tVec, vNeg3Std, 'b');
plot(tVec, vPos3StdSmoothed, 'r--','LineWidth', 0.8);
plot(tVec, vNeg3StdSmoothed, 'r--', 'LineWidth', 0.8);
xlabel("time (s)");
ylabel("v_x error");
lgnd = legend('unsmoothed', 'smoothed', ' \pm 3-\sigma unsmoothed', ...
       '\pm 3-\sigma unsmoothed', '\pm 3-\sigma smoothed', '\pm 3-\sigma smoothed');
lgnd.Location = 'southeast';
lgnd.BoxFace.ColorType='truecoloralpha'
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.55]')
lgnd.NumColumns = 2;
lgnd.Orientation = 'horizontal';
hold off;

sgtitle("Position and Velocity Errors (Unsmoothed vs Smoothed)");

%% Problem 3
% Theory Problem: 
% 1. You are asked to compute the best estimated trajectory at the end of a mission.
%    What algorithm do you select?
% Ans: The best estimated trajectory at the end of a mission will most
%      likely be a batch estimation because it is faster to calculate and
%      is just as accurate as the smoother. The smoother is good if we have
%      a lot of process noise to experience, but as the noise tends to zero
%      it behaves like the batch estimator. 
% 2. Your company has developed a new ranging system and they want you to
%    develop a flight software routine to autonomously compute the spacecraft’s
%    position using that system. What algorithm do you select?
% Ans: For this we would want to use an Extended Kalman Filter and do a
%      post-flight analysis with a smoother. Since this is a new system
%      there will most likely be a lot of process noise so the smoother
%      will handle that better than the Batch Estimation algorithm would.
