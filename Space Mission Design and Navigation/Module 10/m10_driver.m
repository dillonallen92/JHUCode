clear, clc, close all;

% Dillon Allen
% 4/19/2021


% Note: I emailed about some conceptual issues on problem 3/4/5 but I dont
% think I will hear back by the time this is due so I am just going with my
% gut. I explained myself in the email so just let me know what I did wrong
% when you get a chance. Thank you!

% This script is meant to be the driver script requested in the problem set
% for problems 2-5. Added the plot and explanation to problem 1 here as
% comments because it should be easy to fit.

%% Problem 1

r1_IVec = [-1;0];
r2_IVec = [1;0];

% xVec = (0,0)

xVec1 = [0;0];
xVec2 = [1;1];

plot(r1_IVec(1), r1_IVec(2),'b*'); hold on;
plot(r2_IVec(1), r2_IVec(2),'b*'); 
plot(xVec1(1), xVec1(2),'r*');
plot(xVec2(1), xVec2(2),'r*');
plot(0,0,'k*');

% Explanation:
% Since we are doing az/el observations, we can see that xVec1 is
% observable for both I1 and I2. At xVec2, I1 can observe xVec2 but for the
% tangent part of the calculation, we have our angle being 90 degrees above
% I2. This gives an undefined value and therefore we have an issue with
% observability.


%% Problem 2

clear, clc, close all;

% We have two range measuring stations, where one is at (-1,0) and the other is at
% (0,1). Our state is just the 2D position, no velocity.

r1_IVec = [-1;0];
r2_IVec = [1;0];

% xVec = (0,0)

xVec1 = [0;0];
Phi = eye(2,2);

subplot(1,2,1)
circlePlot(xVec1, r1_IVec, r2_IVec);
H1 = findHMat(xVec1, r1_IVec);
Q1 = [H1; H1*Phi];
rank1 = rank(Q1);

% So we have a unique solution where the circles touch, but the rank says
% otherwise

% xVec = [1 1];
xVec2 = [1;1];
subplot(1,2,2)
circlePlot(xVec2, r1_IVec, r2_IVec);
H2 = findHMat(xVec2, r2_IVec);
Q2 = [H2; H2*Phi];
rank2 = rank(Q2);
% The circles intersect at two points so I think we cannot find a unique
% solution. Unfortunately my rank is 1, which confirms it is not observable
% and I am supposed to be getting a contradiction here. I cannot find my
% error though. 

%% Problem 3
clear, clc, close all;

load('HW10_Prob3.mat');

figure;
plot3(rI1_Vec(1), rI1_Vec(2), rI1_Vec(3),'b*'); hold on;
plot3(rI2_Vec(1), rI2_Vec(2), rI2_Vec(3) ,'b*');
plot3(rI3_Vec(1), rI3_Vec(2), rI3_Vec(3) ,'b*');
plot3(xVec(1), xVec(2), xVec(3),'r*');
plot3(0,0,0,'k*');
axis equal;xlabel('x');ylabel('y');zlabel('z'); grid on


opts1.r_IVec = rI1_Vec;
opts2.r_IVec = rI2_Vec;
opts3.r_IVec = rI3_Vec;

[h1, H1, ~] = Provided_measFxnRange_ECI_rStatECI(xVec, w, t, opts1);
[h2, H2, ~] = Provided_measFxnRange_ECI_rStatECI(xVec, w, t, opts2);
[h3, H3, ~] = Provided_measFxnRange_ECI_rStatECI(xVec, w, t, opts3);

H_fullState = [H1; H2; H3];
H_truncated = H_fullState(:,1:3);

Phi_fs = eye(6,6);
Phi_truncated = eye(3,3);

for i = 1:length(xVec)
    if i == 1
        Q_fs = H_fullState;
        Q_truncated = H_truncated;
    else
        Q_fs = [Q_fs; H_fullState*Phi_fs^(i-1)];
        Q_truncated = [Q_truncated; H_truncated * Phi_truncated^(i-1)];
    end
end

rankQfs = rank(Q_fs); % Should be 6, but instead 3. Inst. Pos & Inst. Vel Not observable.
rankQtruncated = rank(Q_truncated); % Should be 3, which is confirmed. Inst. Position is observable.

%% Problem 4
clear, clc, close all;

load('HW10_Prob4And5.mat');

% Plot
figure;
labels = ["station 1", "station 2", "sat", "origin"];
plot3(r1VecHist(1,1), r1VecHist(2,1), r1VecHist(3,1),'b*'); hold on;
text(r1VecHist(1,1), r1VecHist(2,1), r1VecHist(3,1), labels(1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
plot3(r2VecHist(1,1), r2VecHist(2,1), r2VecHist(3,1),'b*');
text(r2VecHist(1,1), r2VecHist(2,1), r2VecHist(3,1), labels(2), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
plot3(xVecHist(1,1), xVecHist(2,1), xVecHist(3,1),'r*');
text(xVecHist(1,1), xVecHist(2,1), xVecHist(3,1), labels(3), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
plot3(0,0,0,'k*');
text(0,0,0,labels(4), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
title("Inertial Placement of Satellite and Range Stations");

% Observability Calculations
H_s1 = zeros(size(xVecHist));
H_s2 = zeros(size(xVecHist));
h_s1 = zeros(size(tspan_sec));
h_s2 = zeros(size(tspan_sec));

for i = 1 : length(tspan_sec)
    xVec = xVecHist(:,i);
    w = 0;
    opts1.r_IVec = r1VecHist(:,i);
    opts2.r_IVec = r2VecHist(:,i);
    [h_s1(:,i), H_s1(:,i), ~] = Provided_measFxnRange_ECI_rStatECI(xVec,w,tspan_sec(i), opts1);
    [h_s2(:,i), H_s2(:,i), ~] = Provided_measFxnRange_ECI_rStatECI(xVec, w, tspan_sec(i), opts2);
end

% Memoize the input
PhiProd = cell(size(Phi));
for i = 1:length(Phi)
   if i == 1
      PhiProd{1,1} = Phi{1,1}; 
   else
       PhiProd{1,i} = Phi{1,i}*PhiProd{1,i-1};
   end
end

Q_s1 = zeros(length(H_s1),length(Phi));
Q_s2 = zeros(length(H_s2),length(Phi));

for i = 1:length(H_s1)
   if i == 1
       Q_s1(i,:) = H_s1(:,i)';
       Q_s2(i,:) = H_s2(:,i)';
   else
       Q_s1(i,:) = H_s1(:,i)'*PhiProd{1,i-1};
       Q_s2(i,:) = H_s2(:,i)'*PhiProd{1,i-1};
   end
end

% Check Ranks
rankQ1 = rank(Q_s1);
if (rankQ1 == length(Phi))
    disp('p4 Station 1 is observable');
else
    disp('p4 Station 1 is not observable');
end

rankQ2 = rank(Q_s2);
if (rankQ2 == length(Phi))
    disp('p4 Station 2 is observable');
else
    disp('p4 Station 2 is not observable');
end

%% Problem 5

% Now calculate the controllability matrix 

Qc = [ Gamma{1,6}, ...
       Phi{1,6}*Gamma{1,5}, ...
       Phi{1,6}*Phi{1,5}*Gamma{1,4}, ...
       Phi{1,6}*Phi{1,5}*Phi{1,4}*Gamma{1,3}, ...
       Phi{1,6}*Phi{1,5}*Phi{1,4}*Phi{1,3}*Gamma{1,2}, ...
       Phi{1,6}*Phi{1,5}*Phi{1,4}*Phi{1,3}*Phi{1,2}*Gamma{1,1}];

rankQc = rank(Qc); % Since the rank is 6, which is the size of our state, it is controllable.

%% Functions used

% Problem 2 functions
function circlePlot(xVec, r1_IVec, r2_IVec)
    dist1 = norm((xVec - r1_IVec));
    dist2 = norm((xVec - r2_IVec));

    % Plot
    plot(r1_IVec(1), r1_IVec(2),'b*'); hold on;
    plot(r2_IVec(1), r2_IVec(2),'b*');
    plot(xVec(1), xVec(2),'r*');
    plot(0,0,'k*');
    viscircles([r1_IVec(1) r1_IVec(2)], dist1, 'Color', 'g');
    viscircles([r2_IVec(1), r2_IVec(2)], dist2, 'Color', 'g');
    axis equal;
    grid on;
end

function H = findHMat(xVec, rI)
    % assuming 2D, so not a general function
    rho = norm((xVec - rI));
    H = [(xVec(1) - rI(1))/rho, (xVec(2) - rI(2))/rho];   
end