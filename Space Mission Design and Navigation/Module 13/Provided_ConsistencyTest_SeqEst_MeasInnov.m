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
% It follows the formulation of Bar-Shalom, Li, and Kirubarajan's
% Estimation with Applications to Tracking and Navigation, Ch 5.4.1-5.4.2
% 
% Core Algorithm:  
%       a normalized innovations test
%                   nu = z - h(xbar)
%                   S = H Pbar H' + R
%                   epsilon_nu = nu' inv(S) nu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) nuVecMat - The [m x k] column vector of measurement
% residuals/innovations for each time k.  It is assumed that the first
% column is a column of zeros.
% 
% 2) SMats - The [m x m x k] matrix of innovation covariance matrices. 
% 
% 3) plotFlag - flag to plot(1) or not (0).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) epsilonVec - the [1 x k] summed and normalized innovations for each
% time-step.
% 
% 2) epsilonBarNu - the sum of epsilonVec, divided by the number of
% measurements (k).
% 
% 3) Exp_epsilon - the expected value of epsilon, assuming the same number of 
% measurements are present for every time-step.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) Sequential estimator.
% 2) Each measurement time produces a full measurement set.
% 3) Stationary statistics
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [epsilonVec, epsilonBarNu, Exp_epsilon] =...
    Provided_ConsistencyTest_SeqEst_MeasInnov(nuVecMat, SMats, plotFlag)

% Sizing:
k = size(nuVecMat,2)-1;

% Sanity check:
if (sum(abs(nuVecMat(:,1)))~=0)
    error('This function expects the first column to be zeros. (nu is stored as k+1, so min storage index is 2)');
end

% Measurement Consistency:
epsilonVec = zeros(1,k);
for ii = 2:k+1
    epsilonVec(ii-1) = nuVecMat(:,ii)' * inv(SMats(:,:,ii)) * nuVecMat(:,ii);
end

epsilonBarNu = (1/k)*sum(epsilonVec);
Exp_epsilon = size(nuVecMat,1);

if (plotFlag==1)
    % Statistics:
    avgDOF = size(nuVecMat,1);
    alpha = 0.025;
    rLower = alpha/2;
    rUpper = 1-(alpha/2);
    boundLower = chi2inv(rLower,avgDOF)
    boundUpper = chi2inv(rUpper,avgDOF)
    %
    % Histogram:
    figure('name','Chi-Square (Meas) Histogram- Bar Chart');
    hist(epsilonVec,linspace(rLower,30*rUpper,100));hold on
    title('\chi^2 (Meas) Histogram, Bar Chart');
    %
    [N,X] = hist(epsilonVec,linspace(rLower,30*rUpper,30));%hist(epsilonVec);
%     subplot(211),bar([boundLower,boundUpper],repmat([0,max(N)],2,1),0.1)
    plot([boundLower,boundLower],[0,max(N)],'r-','linewidth',2)
    plot([boundUpper,boundUpper],[0,max(N)],'r-','linewidth',2)
    %
    figure('name','Chi-Square (Meas) Histogram');plot(X,N/sum(N),'b*')
    title('\chi^2 Test - Measurement')
    Y = chi2pdf(X,size(nuVecMat,1));
    hold on;plot(X,Y,'g-')
    legend('Innovation Data','Expected PDF')
    xlabel('q')
    ylabel('Density, or rough approx from data (not perfect)')
    
    % Bounds:
    figure('name','Chi Square Meas Error Test')
    plot(epsilonVec,'b-*')
    xlabel('Index')
    ylabel('\chi^2 Test')
    hold on
    plot(boundLower*ones(size(epsilonVec)),'r-')
    plot(boundUpper*ones(size(epsilonVec)),'r-')
    legend('\chi^2','Lower 95% bounds (full run)',...
        'Upper 95% bounds (full run)')
    title('\chi^2 Meas Error Test')
    %
    % Check bounds:
    [a,b] = find(epsilonVec>=boundLower);
    ValAboveBound = epsilonVec(b);
    [a,b] = find(ValAboveBound<=boundUpper);
    ValInBounds = ValAboveBound(b);
    numWithinBounds = length(ValInBounds);
    percentWithinBounds = (numWithinBounds/k)*100
end







