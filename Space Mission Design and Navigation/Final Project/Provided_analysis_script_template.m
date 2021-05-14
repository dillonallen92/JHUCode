%% Clean house:
clear all
close all
clc

% Reset random number generator to 1 (repeatability)
seedNumber = 1;
rng(seedNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISSION DESIGN CODE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all;

% set the intitial orbital parameters
figNum=1;                   % figure number for plotting
rE=6378.15;                 % radius of the Earth (km)
muE=398600.4418;            % graviational parameter of the Earth (km^3/s^2)
delT=100;                   % time increment (seconds)
finalTime=7*86400;          % final time ~7 days (seconds)
nB=2*pi/(86400);            % rotation rate of the Earth

% Vary These
a = rE+22000;               % semi-major axis (km)
ecc = 0.7;                  % eccentricity
inc = 85*pi/180;            % inclination (rad) 
w = -90*pi/180;             % argument of peripsis (rad)
Om = 0*pi/180;              % right ascension of the ascending node (rad) 
theta = 0*pi/180;           % true anomaly (rad) 


%ground station locations
% GS(1) - Cordoba, Argentina (COA)
% Keeping this constant
GS(1).lat=-31.416668 * pi/180;
GS(1).long=-64.183334 * pi/180;
statCont.r1 = 'Non-US';

% GS(2) - Anchorage, Alaska 
% Can also be Berlin, Germany

GS(2).lat = 66.160507 * pi/180;
GS(2).long = -153.369141 * pi/180;
statCont.r2 = 'US';

% GS(2).lat = 52.520008 * pi/180;
% GS(2).long = 13.404954 * pi/180;
% statCont.r2 = 'Non-US';


% parameters for ionosphere instrument
alt0=20000;                 % maximize the time above this altitude (km)
lat0=70*pi/180;             % maximize the time above this latitude (rad)


% parameters for South Atlantic Anomaly 
latLb=-10*pi/180;           % maximimize time above this latitude (rad)
latUb=10*pi/180;            % maximize time below this latitude (rad)

longLb=-10*pi/180;          % maximize time above this longitude (rad)
longUb=10*pi/180;           % maximize time below this longitude (rad) 



% convert the orbital elements into a six-state
% Function Completed 5.10.2021
[rV,vV]=orbElToState(a,ecc,inc,Om,w,theta,muE);    %reuse from homework or write this function
state0=[rV;vV];             %initial state

% propagate the s/c for a week
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
tVec=0:delT:finalTime;
[t,x] = ode113(@(t,x)twoBody(t,x,muE),tVec,state0,options);   % use the eoms given with the assignment


% make some plots to visualize the orbit, these funciton will be provided---
index=1000;   % index of the time to plot the position of the spacecraft and ground stations
plotOrb(x,figNum)
plotTimeStamp(x,t,index,rE,nB,GS,figNum)
%-------------------------------------------------------------------------

posIn=x(:,1:3);                                             % convert the inertial position history to the fixed frame
posFixed=posHistInToPosHistFixed(posIn,t,nB);               % use the function provided with the function

%???
[altV,latV,longV]=posFixedToAltLatLong(posFixed,rE);        % write this yourself (page 39)


% These function are given------------------------------------------
% ionisphere parameters
binAltV=aboveValBin(altV,alt0);
binLatV=aboveValBin(latV,lat0);

% times when both contraints are satisfied
binIon=binAltV.*binLatV;

% find the total time when both constraints are satisfied
totTimeIonSat=sum(binIon)*delT;
mdPointsIon=totTimeIonSat/3600


% parameters for South Atlantic Anomaly
binLatSaV=aboveBelowBin(latV,latLb,latUb);
binLongSaV=aboveBelowBin(longV,longLb,longUb);

% times when both contraints are satisfied for South Atlantic Anomaly
binSa=binLatSaV.*binLongSaV;
totTimeSa=sum(binSa)*delT;
mdPointsSa=totTimeSa/3600
%--------------------------------------------------------------------------------


% visually check when contacts exist angV<90 deg.
angV1=getAngStation(GS(1),posFixed);
angV2 = getAngStation(GS(2), posFixed);
figure(2);
hold on;
hp=plot(tVec/86400,angV1*180/pi);
hp.LineWidth=3;
sp = plot(tVec/86400, angV2*180/pi);
sp.LineWidth = 3;
ax=gca;
ax.FontSize=20;
ax.FontName='Times New Roman';
ax.XLabel.String='Time (day)';
ax.YLabel.String='Station Angle (deg)';
ax.XLabel.Interpreter='latex';
ax.YLabel.Interpreter='latex';
box on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NAVIGATION CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%
% Design Knobs:
%%%%%%%%%%%%%%%

%STEP 1:  choose a navSensor statetype
% navSensors.stateType = 'GPS'; measFxnHand=@Provided_measFxnNavSol_ECI;   
navSensors.stateType = 'AWS'; measFxnHand=@Provided_measFxnRange_ECI_multi;  
% navSensors.stateType = 'DSN'; measFxnHand=@Provided_measFxnRange_ECI_multi;
% navSensors.stateType = 'RadioShack'; measFxnHand = @Provided_measFxn_TDOA_and_FDOA_ECEFStat_Multi;
%
%???  Get performance from project assignment for measurement noise
%performance, RMat=???.  

% For AWS
measurementNoiseAWS = 15^2;

% For DSN
measurementNoiseDSN = 1;

% If multiple stations are used, then RMAT needs to grow (diagonally) to
% match.  Write this yourself.

% For two AWS Stations
if isequal(navSensors.stateType,'AWS')
    RMat = measurementNoiseAWS * eye(2);
elseif isequal(navSensors.stateType,'DSN')
    RMat = measurementNoiseDSN * eye(2);
else
    disp('Not a valid station?(ignoring GPS)');
end

% For two DSN Stations
% R = measurementNosieDSN * eye(2);

% Mix:
% R = [measurementNoiseAWS, 0; 0, measurementNoiseDSN;];

%STEP 2:  Define measOpts, contains information for the station positions in the
%measurement models in the Earth fixed frame
%
% The provided multi-station measurement functions call the single station
% measurement functions in a loop and build up h,H,ID.  The input/outputs
% are the same as earlier in the class, with the exception of the
% measurement options, which are for a SPECIFIC station.  Therefore, the
% multi-station functions require a cell array with each element being the
% single station measurement options.  
%   Here is a 2 station example for the range system:
%     measOpts_stat1.r1_ECEFVec  = [1;2;3];
%     measOpts_stat1.tFrameAlign = 0;
%     measOpts_stat2.r1_ECEFVec  = [3;4;5];
%     measOpts_stat2.tFrameAlign = 0;
%     opts_multiStation{1}       = measOpts_stat1;  %define cell and store 1
%     opts_multiStation{2}       = measOpts_stat2;  %store to cell entry 2
%     measOpts                   = opts_multiStation; %to match naming below.

%STEP3:  choose clock

% clockType = "VCTCXO";
% clockType = "OCXO";
clockType = "Rubidium";
% clockType = "Cesium";
%
%???  Get performance from project assignment for clock performance and characteristics



%%%%%%%%%%%%%%
% Conversions:
%%%%%%%%%%%%%%

% transform GS(??) into ECEF coordinates (in units of meters) and store 
% into 'gndStationVecs'. Each column of three represents the position of 
% a single ground station

GS1Lat = GS(1).lat;
GS1Long = GS(1).long;
GS2Lat = GS(2).lat;
GS2Long = GS(2).long;
% Since we are converting ground stations we will use earths radius
rEmeters = rE*1000;
[x1,y1,z1] = lla2ecef(GS1Lat, GS1Long, rEmeters);
gndStationVecs = zeros(3,2);
gndStationVecs(:,1) = [x1;y1;z1];
[x2, y2, z2] = lla2ecef(GS2Lat, GS2Long, rEmeters);
gndStationVecs(:,2) = [x2;y2;z2];
% Nav code uses meters, not kilometers, so convert:
x0_True=zeros(size(x));
x0_True(:,1:3)=x(:,1:3)*1000;    %convert the inertial position to meters
x0_True(:,4:6)=x(:,4:6)*1000;    %convert the inertial position to meters/second
x0_True = x0_True(1,:)';

% measOpts
measOpts_stat1.rStationECEF_m  = gndStationVecs(:,1);
measOpts_stat1.tFrameAlign = 0;
measOpts_stat2.rStationECEF_m = gndStationVecs(:,2);
measOpts_stat2.tFrameAlign = 0;
opts_multiStation{1}       = measOpts_stat1;  %define cell and store 1
opts_multiStation{2}       = measOpts_stat2;  %store to cell entry 2
measOpts                   = opts_multiStation; %to match naming below.


%% Truth-Model Simulation:

% Clock Parameters:
clkSigma0_phase    = 10e-9;
clkSigma0_freq     = 1e-9;  % Ryan updated, 5.12.2021
clkSigma0_freqRate = 1e-15; % effectively 0.
% Scalar values for tuning:
stateScalar_m   = 1e3; % m will be squared below
stateScalar_mps = 1; % m/s will be squared below
procNoiseScalar = 1e-8; % will be squared below
% Number of measurements times:
nk = length(tVec);

% Setup dynamics:
dynFxnHand      = @Provided_dynFxnMod_Earth2Body_ECI;
dynOpts.odeOpts = odeset('RelTol',1e-10,'AbsTol',1e-10);
QMat            = eye(3)*(procNoiseScalar^2);
%
% TM call:
[zVecHist,zNoNoiseVecHist,xTrueVecHist,IDsCells] =...
    Provided_TruthModelSimulationFull_Cells(dynFxnHand,dynOpts,measFxnHand,...
    measOpts,x0_True,tVec,RMat,QMat);


%% Navigation:
nX         = size(x0_True,1);
PBarMat0   = blkdiag( eye(3)*(stateScalar_m^2), eye(3)*(stateScalar_mps^2) );
xStar_0Vec = x0_True + (chol(PBarMat0)*randn(nX,1));
%
YCells = zVecHist;
for ii = 1:nk
    RCells{ii} = RMat;
    QCells{ii} = QMat;
end
% Filter:
[xStarCells,PiCells,nuVecCells,SCells,debug] =...
    Provided_CompAlgForExtSeqProcWSNC_Cells(YCells,IDsCells,RCells,...
    xStar_0Vec,PBarMat0,tVec,dynFxnHand,dynOpts,measFxnHand,measOpts,...
    QCells);


%% Errors:
initialError = norm(xStar_0Vec - xTrueVecHist{1});
finalError   = norm(xStarCells{end} - xTrueVecHist{end});
%
% Consistency Checks:
%
% chi-square measurements test
% ???
% TODO: Ask if okay to use provided Chi-Squared functions...
% function [epsilonVec, epsilonBarNu, Exp_epsilon] =...
%    Provided_ConsistencyTest_SeqEst_MeasInnov(nuVecMat, SMats, plotFlag)
%

plotFlag = 1;
Provided_ConsistencyTest_SeqEst_MeasInnov(nuVecCells, SCells, plotFlag);

% State error:
for ii = 1:nk
    if isempty(xStarCells{ii})
        xeMat(:,ii) = NaN;
    else
        xeMat(:,ii) = xStarCells{ii} - xTrueVecHist{ii};
        PiMats(:,:,ii) = PiCells{ii};
    end
end
% chi-square state error test
% ???
% function [epsilonVec, epsilonBarNu, Exp_epsilon] =...
%    Provided_ConsistencyTest_SeqEst_StateError(xeVecMat, PMats, plotFlag)
Provided_ConsistencyTest_SeqEst_StateError(xeMat, PiMats, plotFlag)

%% Pass times: (visual inspection)
gndStatVisb = zeros(1,nk);
for ii = 1:nk
    if (isempty(nuVecCells{ii}) == 1)
        gndStatVisb(ii) = 0;
    else
        if isequal(measFxnHand,@Provided_measFxn_TDOA_and_FDOA_ECEFStat_Multi)
            gndStatVisb(ii) = numel(nuVecCells{ii})/2;
        elseif isequal(measFxnHand,@Provided_measFxnRange_ECI_multi)
            gndStatVisb(ii) = numel(nuVecCells{ii});            
        end
    end
end
figure('name','For Score - # stations visible & Gnd Contacts')
plot(tVec,gndStatVisb)
hold on;
days_sec = [0:7]*86400;
% Plot a red line each day:
for ii = 1:numel(days_sec)
    line([days_sec(ii),days_sec(ii)], ...
        [0, max(gndStatVisb)],'Color','red','LineStyle','--')
end





%% Position & Velocity Analysis: 

% Compute Sigmas:
for ijk = 1:nk
    SigmasFilt(:,ijk) = diag(PiMats(:,:,ijk).^0.5);
end

%%%%%%%%%%%%%%%%%
% Pos&Vel Errors:
%%%%%%%%%%%%%%%%%
%
nSig = 3;
%
figure('name', 'Position Error');
subplot(411),plot(tVec,xeMat(1,:),'r-')
title('X_{ECI} Position Error')
hold on;
subplot(411),plot(tVec,SigmasFilt(1,:)*nSig,'b-')
subplot(411),plot(tVec,-SigmasFilt(1,:)*nSig,'b-')
xlabel('Time [sec]')
ylabel('Error [m]')
subplot(412),plot(tVec,xeMat(2,:),'r-')
title('X_{ECI} Position Error')
hold on;
subplot(412),plot(tVec,SigmasFilt(2,:)*nSig,'b-')
subplot(412),plot(tVec,-SigmasFilt(2,:)*nSig,'b-')
title('Y_{ECI} Position Error')
xlabel('Time [sec]')
ylabel('Error [m]')
subplot(413),plot(tVec,xeMat(3,:),'r-')
title('X_{ECI} Position Error')
hold on;
subplot(413),plot(tVec,SigmasFilt(3,:)*nSig,'b-')
subplot(413),plot(tVec,-SigmasFilt(3,:)*nSig,'b-')
title('Z_{ECI} Position Error')
legend('Error',['+/- ',num2str(nSig),' \sigma bounds'])
xlabel('Time [sec]')
ylabel('Error [m]')
subplot(414),plot(tVec,gndStatVisb,'k-')

figure('name', 'Velocity Error');
subplot(411),plot(tVec,xeMat(4,:),'r-')
title('X_{ECI} Velocity Error')
hold on;
subplot(411),plot(tVec,SigmasFilt(4,:)*nSig,'b-')
subplot(411),plot(tVec,-SigmasFilt(4,:)*nSig,'b-')
xlabel('Time [sec]')
ylabel('Error [m/s]')
subplot(412),plot(tVec,xeMat(5,:),'r-')
title('Y_{ECI} Velocity Error')
hold on;
subplot(412),plot(tVec,SigmasFilt(5,:)*nSig,'b-')
subplot(412),plot(tVec,-SigmasFilt(5,:)*nSig,'b-')
xlabel('Time [sec]')
ylabel('Error [m/s]')
subplot(413),plot(tVec,xeMat(6,:),'r-')
title('Z_{ECI} Velocity Error')
hold on;
subplot(413),plot(tVec,SigmasFilt(6,:)*nSig,'b-')
subplot(413),plot(tVec,-SigmasFilt(6,:)*nSig,'b-')
legend('Error','+/- 1 \sigma bounds')
xlabel('Time [sec]')
ylabel('Error [m/s]')
subplot(414),plot(tVec,gndStatVisb,'k-')


%%%%%%%%%%%%%%%%
% Position Norm:
%%%%%%%%%%%%%%%%
%
% Compute Sigmas:
xePosNorm_m = ( (xeMat(1,:).^2) + (xeMat(2,:).^2) + (xeMat(3,:).^2) ).^0.5;
xeVelNorm_m = ( (xeMat(4,:).^2) + (xeMat(5,:).^2) + (xeMat(6,:).^2) ).^0.5;
sigmaPosNorm_m = ( (SigmasFilt(1,:).^2) + (SigmasFilt(2,:).^2) + (SigmasFilt(3,:).^2) ).^0.5;
sigmaVelNorm_m = ( (SigmasFilt(4,:).^2) + (SigmasFilt(5,:).^2) + (SigmasFilt(6,:).^2) ).^0.5;
%
figure('name','For Score - Sigmas and Error')
subplot(211);plot(tVec,sigmaPosNorm_m*nSig,'r--'); hold on
subplot(211);plot(tVec,abs(xePosNorm_m),'b-*')
xlabel('Time [sec]')
ylabel('Position Error & \sigma [m]')
title('Filter Position Error and \sigma Plots');
grid on
subplot(212);plot(tVec,sigmaVelNorm_m*nSig,'r--'); hold on
subplot(212);plot(tVec,abs(xeVelNorm_m),'b-*')
legend(['+/- ',num2str(nSig),' \sigma bounds'],'|Pos. Er.|')
xlabel('Time [sec]')
ylabel('Velocity Error & \sigma [m/s]')
title('Filter Velocity Error and \sigma Plots')
grid on

%%%%%%%%%%%
% 3D orbit:
%%%%%%%%%%%
for ii = 1:nk
   xVecTrue(:,ii) = xTrueVecHist{ii};
   xVecEst(:,ii)  = xStarCells{ii};
end
figure('name','3D Orbit, True & Est')
plot3(xVecTrue(1,:),xVecTrue(2,:),xVecTrue(3,:),'g-*')
hold on; grid on; axis equal
plot3(xVecEst(1,:),xVecEst(2,:),xVecEst(3,:),'r-*')
xlabel('x')
ylabel('y')
zlabel('z')

%% Clocks Analysis:

% Simulate the dynamics of your clock and compare to requirements.
%   You really only care about the covariance propagation, not the actual
%   simulation of the phase/frequency/frequency rate time-histories.  
%
% ???
%
% Check if IDCells are empty... If so, propagate the covariance. Otherwise,
% reset until empty again
% Initial State Vector
% clkSigma0_phase    = 10e-9;
% clkSigma0_freq     = 1e-9;  % Ryan updated, 5.12.2021
% clkSigma0_freqRate = 1e-15; % effectively 0.
% P covariance Value: phase_error = 10ns, freq_error = 1e-9 again
x0Clock = [clkSigma0_phase; clkSigma0_freq; clkSigma0_freqRate];
PClockInit = diag(x0Clock);
PClockCells{1} = PClockInit;
uClock = 0; vClock = 0;

% Clock QMat
switch(clockType)
    case 'VCTCXO'
        sigma1clock = 2.24e-11;
        sigma2clock = 4.44e-10;
        sigma3clock = 0;
    case 'OCXO'
        sigma1clock = 1.12e-11;
        sigma2clock = 7.04e-11;
        sigma3clock = 0;
    case 'Rubidium'
        sigma1clock = 2.24e-12;
        sigma2clock = 5.06e-13;
        sigma3clock = 0;
    case 'Cesium'
        sigma1clock = 1e-10;
        sigma2clock = 2.81e-14;
        sigma3clock = 0;
    otherwise
        disp('Clock not found?');
end

QMatClock = diag([sigma1clock; sigma2clock; sigma3clock]);
xClockOut = x0Clock;
for idCount = 1 : length(IDsCells)-1
    % dynamics for each time step
    tVecClock = [tVec(idCount), tVec(idCount+1)];
    [tClockOut, xClockOut, PhiClock, GuClock, GvClock] = ... 
        Provided_dynFxnHand_singleClockModel_ZT(xClockOut, uClock, vClock, tVecClock);
    % Check the current IDsCell. If empty, that means no ground contact.
    % If no ground contact, propagate the state from that time to another
    % PClockCells{idCount+1} = PhiClock * PClockCells{idCount} *
    % transpose(PhiClock) + GuClock * QMatClock * transpose(GuClock);
    % else, PClockCells{idCount+1} = PClockInit
    
    if isempty(IDsCells{idCount})
        PClockCells{idCount+1} = PhiClock * PClockCells{idCount} * transpose(PhiClock) + ...
                                 GuClock * QMatClock * transpose(GuClock);
    else
        PClockCells{idCount+1} = PClockInit;
    end    
end

nPClockCells = length(PClockCells);
for ii = 1:nPClockCells
   phaseErrorVec(ii) = sqrt(PClockCells{ii}(1,1));
   freqErrorVec(ii) = sqrt(PClockCells{ii}(2,2));
   freqPhaseErrorVec(ii) = sqrt(PClockCells{ii}(3,3));
end

clockLenVec = 1:nPClockCells;
clockFigTitle = sprintf('%s clock error',clockType);
figure('name', clockFigTitle)
subplot(311)
plot(clockLenVec, phaseErrorVec);
xlabel('index');
ylabel('error');
title('Phase Error over Time');

subplot(312)
plot(clockLenVec, freqErrorVec);
xlabel('index');
ylabel('error');
title('Frequency Error over Time');

subplot(313)
plot(clockLenVec, freqPhaseErrorVec);
xlabel('index');
ylabel('error');
title('Frequency Phase Error over Time');

%% Plotting the solution: (Visual verification that the location are appropriate)
figure('name','2D Earth')
theta = [-179:180];
phi   = [-89:90];
obs   = 1 % 360x180 2d array of data ???
[lat_grid,lon_grid] = meshgrid(theta,phi); % should both be 360x180 to match data
load coast % loads lat and long variables that define the coastline
ax = worldmap('World') % also try axesm as it gives more options
geoshow(lat,long) % draw the coastlines
%
% Plot Ground Stations:
M           = size(gndStationVecs,2);
xECEFMat_km = gndStationVecs'/1e3;
for ii = 1:M
    latitude_deg(ii)  = (180/pi)*atan2(gndStationVecs(3,ii),norm(gndStationVecs(1:2,ii)));
    longitude_deg(ii) = (180/pi)*atan2(gndStationVecs(2,ii),gndStationVecs(1,ii));
end
%
rad = 2;
for ii = 1:M
    [unitVecMat] = computeUnitVectorsOnCircle(20, 2);
    Lats = latitude_deg(ii) + (unitVecMat(1,:)*rad);
    Lons = longitude_deg(ii) + (unitVecMat(2,:)*rad);  
    % Plot:
    geoshow(ax, Lats, Lons,...
        'DisplayType', 'polygon', 'FaceColor', [.45 .60 .30])
end
figure('name','XYZ plot for stations')
plot3(xECEFMat_km(:,1),xECEFMat_km(:,2),xECEFMat_km(:,3),'g*')
axis equal; hold on; grid on
plot3(0,0,0,'r*')
    
    
%% Cost of solution:
% Ground Station Cost (per)
costAWS = 100000;
costDSN = 400000;
costHumanOp = 100000;

% Comm Cost
costCommUS = 10000;
costCommNon_US = 100000;
costCommAntartica = 500000;

% Clock Costs
costVCTCXO = 0;
costOCXO = 10000;
costRubidium = 500000;
costCesium = 1000000;

% Initialize Cost
cost = 0;
if isequal(navSensors.stateType,'AWS')
    cost = cost + costHumanOp + costAWS*2;
end

if isequal(navSensors.stateType,'DSN')
    cost = cost + costHumanOp + costDSN*2;
end

if isequal(statCont.r1,'Non-US')
   cost = cost +  100000;
else
   cost = cost + 10000;
end

if isequal(statCont.r2,'Non-US')
       cost = cost +  100000;
else
   cost = cost + 10000;
end

switch(clockType)
    case 'VCTCXO'
        cost = cost + costVCTCXO;
    case 'OCXO'
        cost = cost + costOCXO;
    case 'Rubidium'
        cost = cost + costRubidium;
    case 'Cesium'
        cost = cost + costCesium;
    otherwise
        disp('Clock not found?');
end

fprintf("Total cost of the system: %d\n", cost);