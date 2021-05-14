%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ryan Mitch
% Date Modified: 2-2-2019
% Copyright (c) 2019 Ryan Mitch and The Johns Hopkins University.  All
% rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same algorithm as CompAlgForExtSeqProcWSNC_withNoiseLogic.m, but
% trimmed for distribution.
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
%  [h,H,id] = measFxnHandTest3(x,n,t,opts)
%   where x is the state after propagation, n is the noise vector (should
%   be 0 in estimation, but non-zero in truth-model simulation, t is the
%   measuremnt time, and opts are the measurement options.
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
% 1) xStarCells - The cell array of the estimated states resulting
% from the extended sequential estimator with state noise compensation.
%
% 2) PMat - The cell array of the matrix of the state estimate covariances.  
% 
% 3) nuVecMat - The cell array of measurement innovations. 
%  
% 4) SMats - The cell array of covariances mapped to the measurement space.
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
function [xStarCells,PiCells,nuVecCells,SCells,debug] =...
    Provided_CompAlgForExtSeqProcWSNC_Cells(YCells,IDsCells,RCells,...
    xStar_0Vec,PBarMat0,tVec,dynFxnHand,dynFxnOpts,measFxnHand,measFxnOpts,...
    QCells)

%% Optional param: (trimmed down from nominal)
%
debug = [];

%% Filter:
% BLOCK 1:  Initialize at t0
xStar_tim1Vec = xStar_0Vec(:,1); 
P_im1         = PBarMat0;
k             = length(tVec)-1;
tim1          = tVec(1); % t0
% Storage:
xStarCells    = {};
xStarCells{1} = xStar_tim1Vec;
PiCells{1}    = P_im1;
%
% Loop over k measurements
for ii = 1:k
    %
    % BLOCK 2:  (A) Read the next observation
    ti   = tVec(ii+1); %[t0,t1,...,tk]
    Yi   = YCells{ii}; 
    Qi   = QCells{ii};
    Ri   = RCells{ii};
    IDsi = IDsCells{ii};
    %
    % BLOCK 3: Integrate reference trajectory and state transition matrix
    % from t_im1 to t_i:
    uk = [];
    [~,xStar_tiVec,Phi_titim1,Gu_titim1,Gv_titim1] = feval(dynFxnHand,...
        xStar_tim1Vec,[],uk,[tim1,ti],dynFxnOpts);
    %
    % BLOCK 4: Time Update
    Pbar_i = (Phi_titim1 * P_im1 * Phi_titim1') + ( Gu_titim1 * Qi * Gu_titim1' );      
    %
    % BLOCK 5: Compute: observation, observation-state matrixu, gain matrix
    [GiVec, HTildeMat, IDs] = feval(measFxnHand,xStar_tiVec,0,ti,measFxnOpts);
    if ~all(size(IDs)==size(GiVec))
        error('Need an ID for every measurement.')
    end
    
    % If not measurements, just continue:
    if isempty(Yi) || isempty(GiVec)
        xStarCells{ii+1} = xStar_tiVec;
        Pi               = Pbar_i;
        xStarCells{ii+1} = xStar_tiVec;
        PiCells{ii+1}    = Pi;
        nuVecCells{ii+1} = [];
        SCells{ii+1}     = [];
        continue
    end
    %
    % Numerically compute Jacobian.  Uncomment if desired (and not done in
    % the measurement model function)
%     [HTildeMat] = computeJacobian_measMod(measFxnHand,xStar_tiVec,[],ti,measFxnOpts);
    %
    % Match:
    %??? match IDs and trim to the common subset
    % Yi        = Yi(???);
    % GiVec     = GiVec(???);
    % HTildeMat = HTildeMat(???,:);
    % Ri        = Ri(???,???);
    
    [IDcommoni,IMeas,IEst] = intersect(IDsi,IDs);   

    Yi        = Yi(IMeas);   

    GiVec     = GiVec(IEst);   

    HTildeMat = HTildeMat(IEst,:);   

    Ri        = Ri(IMeas,IMeas);
    %
    yiVec = Yi - GiVec;
    Si    = (HTildeMat * Pbar_i * HTildeMat') + Ri;
    Ki    = (Pbar_i * HTildeMat') * inv( Si );
    %
    % BLOCK 6: Measurement and Reference Update
    xhat_i = Ki * yiVec;
    L      = eye(size(xStar_tiVec,1)) - Ki * HTildeMat;
    Pi     = (L*Pbar_i*L') + (Ki * Ri * Ki');    
    xStar_tiVec = xStar_tiVec + xhat_i;
    % Update the variable names for next measurement set:
    tim1          = ti;
    xStar_tim1Vec = xStar_tiVec;
    P_im1         = Pi;
    %
    % Store for output:
    xStarCells{ii+1} = xStar_tiVec;
    PiCells{ii+1}    = Pi;
    nuVecCells{ii+1} = yiVec;
    SCells{ii+1}     = Si; 
end
return