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