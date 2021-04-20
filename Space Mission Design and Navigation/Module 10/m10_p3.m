clear, clc, close all;

load('HW10_Prob3.mat');

%% 
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