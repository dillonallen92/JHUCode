clear, clc, close all;
% Problem 1
% Create a delta-V plot vs mass fraction
% for Isp of cold gas (75), monoprop (220)
% biprop (316) and hall effect thrusters (1380)

%Isp
Isp_cg = 75; Isp_mp = 220; Isp_bp = 316; Isp_het = 1380;

% gravitational constant
g0 = 9.81;

% Mass of spacecraft
Msc = 120;

% Remaining fraction
% Mpl + Mpr = 60
% Mpr = 60 - Mpl
% Max payload can be 60kg (0 pr)

Mpl = 0:0.01:60;
Mpr = 60 - Mpl;
Mi = Msc + Mpl + Mpr;
Mf = Msc + Mpl;
y = @(Isp) Isp*g0*log(Mi./Mf);
ycg = y(Isp_cg);
ymp = y(Isp_mp);
ybp = y(Isp_bp);
yhet = y(Isp_het);
x = Mi./Mf;
plot(Mpl,ycg)
hold on
plot(Mpl,ymp);
plot(Mpl, ybp);
plot(Mpl,yhet);
legend('Cold Gas', 'Monoprop', 'Biprop', 'Hall Effect Thursters');
xlabel('Payload Mass (Kg)');
ylabel('\DeltaV (m/s)');
title('\DeltaV vs Payload Mass (Kg)');


