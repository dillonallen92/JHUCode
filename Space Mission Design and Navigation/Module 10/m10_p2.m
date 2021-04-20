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