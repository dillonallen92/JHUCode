clear, clc, close all;

% This script is meant to be the driver script requested in the problem set
% for problems 2-5

%% Problem 2

%% Problem 3

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