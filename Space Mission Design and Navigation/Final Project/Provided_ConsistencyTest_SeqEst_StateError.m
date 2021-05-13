%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 12-1-2019
% Copyright (c) 2019 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comprises the core computations of a Sequential Estimator's
% measurement innovation test. 
% It follows the formulation of ???
% 
% Core Algorithm:  
%       a normalized innovations test
%                   xe = xTrue - xEst
%                   P = posterior covariance
%                   epsilon_xe = xe' inv(P) xe
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) xeVecMat - The [m x k] column vector of state errors for each time k.
% It is assumed that the first column is the initial state error before a
% measurement update occurs. (e.g., x0)
% 
% 2) PMats - The [m x m x k] matrix of state covariance matrices. 
% 
% 3) plotFlag - flag ot plot(1) or not (0).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) epsilonVec - the [1 x k] summed and normalized state errors for each
% time-step.
% 
% 2) epsilonBarNu - the sum of epsilonVec, divided by the number of
% measurements (k). ???
% 
% 3) Exp_epsilon - the expected value of epsilon, assuming the same number 
% of measurements are present for every time-step.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) Sequential estimator.
% 2) Fixed state size.
% 3) Stationary statistics.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [epsilonVec, epsilonBarNu, Exp_epsilon] =...
    Provided_ConsistencyTest_SeqEst_StateError(xeVecMat, PMats, plotFlag)

% Sizing:
k = size(xeVecMat,2)-1;

% Measurement Consistency:
epsilonVec = zeros(1,k);
for ii = 1:k+1
    epsilonVec(ii) = xeVecMat(:,ii)' * inv(PMats(:,:,ii)) * xeVecMat(:,ii);
end

epsilonBarNu = (1/k)*sum(epsilonVec);
Exp_epsilon = size(xeVecMat,1);

if (plotFlag==1)
    [N,X] = hist(epsilonVec);
    figure;plot(X,N/sum(N),'b*')
    Y = chi2pdf(X,size(xeVecMat,1));
    hold on;plot(X,Y,'g-')
    legend('State Error Data','Expected PDF')
    xlabel('q')
    ylabel('Density, or rough approx from data (not perfect)')
    title('\chi^2 Test - State Error')
    
    figure('name','Chi Square State Error Test')
    plot(epsilonVec,'b-*')
    xlabel('Index')
    ylabel('\chi^2 Test')
    %
    avgDOF = size(xeVecMat,1);
    alpha = 0.05;
    rLower = alpha/2;
    rUpper = 1-(alpha/2);
    boundLower = chi2inv(rLower,avgDOF)
    boundUpper = chi2inv(rUpper,avgDOF)
    hold on
    plot(boundLower*ones(size(epsilonVec)),'r-')
    plot(boundUpper*ones(size(epsilonVec)),'r-')
    legend('\chi^2','Lower 95% bounds (full run)',...
        'Upper 95% bounds (full run)')
    title('\chi^2 State Error Test')
    %
    % Check bounds:
    [a,b]         = find(epsilonVec>=boundLower);
    ValAboveBound = epsilonVec(b);
    [a,b]         = find(ValAboveBound<=boundUpper);
    ValInBounds   = ValAboveBound(b);
    numWithinBounds     = length(ValInBounds);
    percentWithinBounds = (numWithinBounds/k)*100
    
end






