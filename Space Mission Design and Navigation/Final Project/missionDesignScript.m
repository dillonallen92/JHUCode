%% Clean house:
clear
close all
clc

% Reset random number generator to 1 (repeatability)
seedNumber = 1;
rng(seedNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISSION DESIGN CODE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the intitial orbital parameters
% figNum=1;                   % figure number for plotting
rE=6378.15;                 % radius of the Earth (km)
muE=398600.4418;            % graviational parameter of the Earth (km^3/s^2)
delT=100;                   % time increment (seconds)
finalTime=7*86400;          % final time ~7 days (seconds)
nB=2*pi/(86400);            % rotation rate of the Earth

mdPointsIont = 0;
mdPointsSat = 0;
count = 1;

for r = 22000 : 1000 : 30000
    for e = 0.2: 0.05: 0.8
        for i = 71:1:120
            for omega = -180 : 10 : 180
                for Omega = -180 : 10 : 180
                    for Theta = -180 : 10 : 180
                
                        a = rE+r;                  % semi-major axis (km)
                        ecc = e;                   % eccentricity
                        inc = i*pi/180;            % inclination (rad) 
                        w = omega*pi/180;          % argument of peripsis (rad)
                        Om = Omega*pi/180;             % right ascension of the ascending node (rad) 
                        theta = Theta*pi/180;          % true anomaly (rad) 


                        %ground station locations
                        % GS(1) - Cordoba, Argentina (COA)
                        GS(1).lat=-31.416668 * pi/180;
                        GS(1).long=-64.183334 * pi/180;

                        % GS(2) - Anchorage, Alaska 
                        GS(2).lat = 66.160507 * pi/180;
                        GS(2).long = -153.369141 * pi/180;


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


                        % % make some plots to visualize the orbit, these funciton will be provided---
                        % index=1000;   % index of the time to plot the position of the spacecraft and ground stations
                        % plotOrb(x,figNum)
                        % plotTimeStamp(x,t,index,rE,nB,GS,figNum)
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
                        mdPointsIon=totTimeIonSat/3600;


                        % parameters for South Atlantic Anomaly
                        binLatSaV=aboveBelowBin(latV,latLb,latUb);
                        binLongSaV=aboveBelowBin(longV,longLb,longUb);

                        % times when both contraints are satisfied for South Atlantic Anomaly
                        binSa=binLatSaV.*binLongSaV;
                        totTimeSa=sum(binSa)*delT;
                        mdPointsSa=totTimeSa/3600;
                        %--------------------------------------------------------------------------------


                        % % visually check when contacts exist angV<90 deg.
                        % angV1=getAngStation(GS(1),posFixed);
                        % angV2 = getAngStation(GS(2), posFixed);
                        % figure(2);
                        % hold on;
                        % hp=plot(tVec/86400,angV1*180/pi);
                        % hp.LineWidth=3;
                        % sp = plot(tVec/86400, angV2*180/pi);
                        % sp.LineWidth = 3;
                        % ax=gca;
                        % ax.FontSize=20;
                        % ax.FontName='Times New Roman';
                        % ax.XLabel.String='Time (day)';
                        % ax.YLabel.String='Station Angle (deg)';
                        % ax.XLabel.Interpreter='latex';
                        % ax.YLabel.Interpreter='latex';
                        % box on;
                        if (mdPointsIon > 0 && mdPointsSa > 0)
                            SolCell(count, :) = [mdPointsIon, mdPointsSa, r, ecc, i, omega, Omega, Theta];
                            count = count + 1;
                        end
                    end
                end
            end
        end
    end
end
                
fprintf('Max Ion: %d, Max SAA: %d\n', max(SolCell(:,1)), max(SolCell(:,2)));
idxSaOver6 = find(SolCell(:,2) > 6);
SolCellSaOver6 = SolCell(idxSaOver6,:);
idxIon30Sa3 = find(SolCell(:,1) > 30 & SolCell(:,2) > 3);
SolCellIon30Sa3 = SolCell(idxIon30Sa3,:);
