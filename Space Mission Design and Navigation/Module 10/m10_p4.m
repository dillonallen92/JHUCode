clear, clc, close all;

load('HW10_Prob4And5.mat');

H_s1 = [];
H_s2 = [];
h_s1 = [];
h_s2 = [];
Q_s1 = [];
Q_s2 = [];
for i = 1 : length(tspan_sec)
    xVec = xVecHist(:,i);
    w = 0;
    opts1.r_IVec = r1VecHist(:,i);
    opts2.r_IVec = r2VecHist(:,i);
    [h_s1(:,i), H_s1(:,i), id] = Provided_measFxnRange_ECI_rStatECI(xVec,w,tspan_sec(i), opts1);
    [h_s2(:,i), H_s2(:,i), id] = Provided_measFxnRange_ECI_rStatECI(xVec, w, tspan_sec(i), opts2);
end

