%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 2-2-2019
% Copyright (c) 2019 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comprises the core computations of a Sequential Estimator.
% This implementation follows flow diagrams of Tapley, Schutz, and Born 
% Chapter 4.7 (2004 edition).
% 
% Core Algorithm:  
%   an Extended Sequential Linear Weighted Least Squares Estimator with 
%   State Noise Compensation.
%   This is a single iteration.
%                   xb = Phi x
%                   Pb = (Phi P Phi') + ( G Q G')
%                   xh = xb + K(y-Hxb)
%                   P = (I - KH)*Pb
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) YMat - The [m x k] column vector of measurements.
% 
% 2) RMats - The [m x m x k] matrix of measurement noise covariance matrix. 
% (this is the inverse of the weighting matrix in WLS)
% 
% 3) xStar_0Vec - The [n x 1] vector of the reference state vector. 
% 
% 4) PBarMat0 - The [n x n] matrix of the prior state estimate covariances.
% (this is the inverse of the WBar matrix in WLS w/ prior).
% 
% 5) tVec - The [(k+1) x 1] vector of times to solve the state for.  The
% first time is the time of the initial state (tm1) and does not correspond
% to a measurement time.  All other times corrspond to a measurement time.
% 
% 6) dynFxnHand - the dynamics model function handle.  This is meant to
% generalize the function.  A standard MATLAB function can be called by
% definining a variable for that function:  
%   functionHandleVariable = @functionName
% In this case the function is supposed to follow this format:
%   [tout,xOut,F,Gu,Gv] = dynFxnHandTest3(x,u,v,tV,opts) 
%   where x is the state vector, u is the disturbance input, v is a
%   control input (0 here), tVec is the time vector [tfrom, tto], and
%   opt is an options function.  
%   tout is the output time (should be tV(end)), xOut is the state after
%   propagation, F is the integrated state transistion matrix, Gu is the
%   process noise influence matrix, Gv is the control influence
%   matrix (0 here).
%  
% 7) measFxnHand - the measurement model function handle.  Similar to
% dynFxnHand, but for the measurements.  In this case it should follow this
% format:
%  [h,H,id] = measFxnHandTest3(x,n,opts)
%   where x is the state after propagation, n is the noise vector (should
%   be 0 in estimation, but non-zero in truth-model simulation.
% 
% 8) QMats - the [n x n x k] process noise covariance matrices for each of
% the k times steps.  The process noise (u) is assumed to be 0 mean (E[u]=0) 
% and have Q covariance (E[(u-ubar)(u-ubar')] = Q).  This needs to be 
% selected at the same time as the dynamcis model, as the dynamics model's
% Gu will propagate that process noise.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) xHatVec - The [n x 1] column vector of the estimated states resulting
% from the extended sequential estimator with state noise compensation.
%
% 2) PMat - The [n x n] matrix of the state estimate covariances.  
% 
% 3) nuVecMat - The [nz x k] measurement innovations. 
%  
% 4) SMats - The [nz x nz x k] covariances mapped to the measurement space.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 1) Linearized system (no iteration or optimization routine) based upon a
% provided initial reference state (xStar0).  This algorithm updates its
% reference at each time step.
% 2) Sequential estimator.
% 3) Each measurement time produces a full measurement set.
% 4) Stationary statistics
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation here follows the TSB book and it is meant to be
% understandable to the reader (for class).  There are more efficient and
% stable implementations of this algorithm.  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xStarMat,PiMats,nuVecMat,SMats,debug] = Provided_CompAlgForExtSeqProcWSNC(YMat,...
    IDsCells,RMats,xStar_0Vec,PBarMat0,tVec,dynFxnHand,dynFxnOpts,...
    measFxnHand,measFxnOpts,QMats, varargin)

%% Optional param:
%
% Needed for some examples (observability).
if ( isempty(varargin) == 1)
    returnHMatsFlag = 0;
    returnPhiMatsFlag = 0;
    returnPiBarMatsFlag = 0;
    returnxBarFlag = 0;
else
    returnHMatsFlag = varargin{1}.returnHMatsFlag;
    returnPhiMatsFlag = varargin{1}.returnPhiMatsFlag;
    returnPiBarMatsFlag = varargin{1}.returnPiBarMatsFlag;
    returnxBarFlag = varargin{1}.returnxBarFlag;
end

debug = [];

%% Filter:
% BLOCK 1:  Initialize at t0
n             = size(xStar_0Vec,1);
xStar_tim1Vec = xStar_0Vec(:,1); 
P_im1         = PBarMat0;
k             = length(tVec)-1;
tim1          = tVec(1); % t0
% Storage:
xStarMat      = zeros(n,k);
PiMats        = zeros(n,n,k);
xStarMat(:,1) = xStar_tim1Vec;
PiMats(:,:,1) = P_im1;
if (returnPiBarMatsFlag == 1)
    debug.Pbar_i{1} = PBarMat0;      
end
%
% Loop over k measurements
for ii = 1:k
    %
    % BLOCK 2:  (A) Read the next observation
    ti = tVec(ii+1); %[t0,t1,...,tk]
    Yi = YMat(:,ii); 
    Ri = RMats(:,:,ii);
    Qi = QMats(:,:,ii);
    %
    % BLOCK 3: Integrate reference trajectory and state transition matrix
    % from t_im1 to t_i:
    [~,xStar_tiVec,Phi_titim1,Gu_titim1,~] = ...
        feval(dynFxnHand,xStar_tim1Vec,[],[],[tim1,ti],dynFxnOpts);
    % handle, x,u,v,tspan,options
    if (returnxBarFlag == 1)
       debug.xBar{ii+1} = xStar_tiVec;        
    end
    %
    % BLOCK 4: Time Update
    Pbar_i = (Phi_titim1 * P_im1 * Phi_titim1') + ...
        ( Gu_titim1 * Qi * Gu_titim1' );        
    %
    % BLOCK 5: Compute: observation, observation-state matrixu, gain matrix
    [GiVec, HTildeiMat, ~] = feval(measFxnHand,xStar_tiVec,0,ti,measFxnOpts);
    % Use numerical derivatives at the moment (slow):
    [HTildeiMat] = computeJacobian_measMod(measFxnHand,xStar_tiVec,[],ti,measFxnOpts);   
    yiVec = Yi - GiVec;
    Si    = (HTildeiMat * Pbar_i * HTildeiMat') + Ri;
    Ki    = (Pbar_i * HTildeiMat') * inv( Si );
    %
    % BLOCK 6: Measurement and Reference Update
    xhat_i        = Ki * yiVec;
    Pi            = (eye(n) - (Ki*HTildeiMat)) * Pbar_i;
    xStar_tiVec   = xStar_tiVec + xhat_i;
    % Update the variable names for next measurement set:
    tim1          = ti;
    xStar_tim1Vec = xStar_tiVec;
    P_im1         = Pi;
    %
    % Store for output:
    xStarMat(:,ii+1) = xStar_tiVec;
    PiMats(:,:,ii+1) = Pi;
    % Allocate space if necessary for storage...
    if (ii==1)
        nuVecMat = zeros(size(yiVec,1),k);
        SMats    = zeros(size(Si,1),size(Si,2),k);
    end
    nuVecMat(:,ii+1) = yiVec;
    SMats(:,:,ii+1)  = Si; 
    %
    % Optional output:
    if (returnHMatsFlag == 1)
       debug.HMats{ii+1} = HTildeiMat;
    end
    if (returnPhiMatsFlag == 1)
       debug.PhiMats{ii+1} = Phi_titim1;        
    end
    if (returnPiBarMatsFlag == 1)
       debug.Pbar_i{ii+1} = Pbar_i;        
    end
end
return