clear, clc, close all;

% Implement a range measurement model but include the effect of light time.
% Use the template to make a function script

%% Part 1
load('HW8_Prob4.mat');

x = [rSat; rdotSat];
opts.derFlag = 0;
opts.tx = x_TxECEF;
t = 0;
w = 0;
zOrh = TwoWayRange(x,w,t,opts);
rho_ideal = zOrh;

%% Light Time edition

c = 3e8; % m/s
error = 1;
itCount = 0;
while 1
    t_a = t + zOrh/c;
    delT = t_a - t;
    r_new = rSat + (delT)*rdotSat;
    x_new = [r_new; rdotSat];
    zOrh_new = TwoWayRange(x_new, w, delT, opts);
    error = zOrh_new - zOrh;
    if error < 1e-6
        break
    else
        zOrh = zOrh_new;
        t = delT;
        itCount = itCount + 1;
    end
end

%% 
