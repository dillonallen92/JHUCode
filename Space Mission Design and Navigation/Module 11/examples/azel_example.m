clear, clc, close all;

% Attempt on the Az/El example as seen on slides 37 - 41
load('NLWLS_Ex.mat');

% Process
% Convert from ECEF to ENV Coordinates
% Compute az/el
% Find Jacobian through numerical differencing
% x*_0 = x_0
% x_hat(1) = x_0 - inv(H'(x_0)*H(x_0))*(H'(x_0)*(h(x_0) - z(1)))
% x_star_i+1 = x_hat(i)
% repeat until norm(x_star(i+1) - x_star(i)) < 1e-6

% Questions:
% where do I get the measurements for the example problem
%   - NLWLS_Ex has the zVec measurement vector
% How do I do this problem with 3 instrument locations
%   - run the measurement model 3 times, make a full [h] matrix and [H]
%     matrix
% The flow chart in 4.6 uses the STM, but the notes have a process without it, which one?
%   - no need for STM because all at single time epoch

r1_IVec = rStation1_m;
r2_IVec = rStation2_m;
r3_IVec = rStation3_m;
rVec = rSat_ECEF_m;
xVec1 = rSat_ECEF_m - r1_IVec;
xVec2 = rSat_ECEF_m - r2_IVec;
xVec3 = rSat_ECEF_m - r3_IVec;
x_star = x_0;
t = 0;
opts.t = 0;
[h1,H1] = AzElMeasurementModel(xVec1,t,opts);
[h2,H2] = AzElMeasurementModel(xVec2,t,opts);
[h3,H3] = AzElMeasurementModel(xVec3,t,opts);
h = [h1; h2; h3];
H = [H1; H2; H3];
zVec = zVec .* 1e6;
x_hat = x_star - inv(H'*WMat * H)*(H' * WMat * (h - zVec ));
x_star = x_hat;
while 1
    [h1,H1] = AzElMeasurementModel(x_star,t,opts);
    [h2,H2] = AzElMeasurementModel(x_star,t,opts);
    [h3,H3] = AzElMeasurementModel(x_star,t,opts);
    h = [h1; h2; h3];
    H = [H1; H2; H3];
    x_hat = x_star - inv(H'*WMat * H)*(H' * WMat * (h - zVec));
    if (norm(x_hat - x_star) < 1e-6)
        disp('converges');
        disp('x_hat');
        x_hat
        break;
    elseif (j > 1000)
            disp("run time exceeded");
            break;
    else
        x_star = x_hat;
        j = j + 1;
    end
end


function [z,H] = AzElMeasurementModel(x,t,opts)

    % Constant
    omegaEarth = 7.2921159e-5;
    deltaT = t - opts.t;
    
    % Measurement equations here:
    % zOrh = ???
    
    x = [x(1) x(2) x(3)]';
    rtVec = RotTrop(x) * (RotZ(omegaEarth, deltaT) * x);
    rt = sqrt(rtVec(1)^2 + rtVec(2)^2 + rtVec(3)^2);
    z = [atan2(rtVec(2), rtVec(1)); asin(rtVec(3)/rt)];
    
    % Numerical Jacobian
    epsilon = 1e-7;
    for j = 1:2
        for i = 1:3
        x(i) = x(i) + epsilon;
        rtVecptrb = RotTrop(x) * (RotZ(omegaEarth, deltaT) * x);
        rtptrb = sqrt(rtVecptrb(1)^2 + rtVecptrb(2)^2 + rtVecptrb(3)^2);
        h(:,i) = [atan2(rtVecptrb(2), rtVecptrb(1)); asin(rtVecptrb(3)/rtptrb)];
        H(:,i) = (h(i) - z)/epsilon;
        x(i) = x(i) - epsilon;
        end
    end
    
    % Helper functions for transforming frames
    function zRotationMat = RotZ(omegaEarth, deltaT)
        theta = omegaEarth * deltaT;
        zRotationMat = [cos(theta), sin(theta), 0;
        -sin(theta), cos(theta), 0;
        0, 0, 1;];
    end
    
    function tropRotMat = RotTrop(x)
        xy = sqrt(x(1)^2 + x(2)^2);
        lambda = acos(x(1)/xy);
        phi = atan2(x(3), xy);
        tropRotMat = [-sin(lambda), cos(lambda), 0;
        -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi);
        cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi)];
    end

end

