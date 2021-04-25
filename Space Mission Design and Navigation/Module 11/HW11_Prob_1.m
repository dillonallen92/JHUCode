clear, clc, close all;

% Problem 1
% Find:
% 1. Initial State Error
% 2. Estimated State from the Least Squares Method
% 3. Final State Error
% 4. Discussion on improvement

load('HW11_Prob_1_data.mat');

initialStateError = norm(xVecTrue - xInitial);
H = [1 0; 1 1; 1 2; 1 4; 1 5];
x_hat = inv(H'*H)*H'*zVec;
finalStateError = norm(xVecTrue - x_hat);
