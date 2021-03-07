clear, clc, close all;

% Consider a 120 degree 3-dB beamwidth (full angle) antenna on a spacecraft
% at 700 km altitude, body mounted to the spacectaft with its  boresight
% pointing nadir.

%% Part A
% Use the antenna gain approximation to estimate the decrease in uplink due
% to beam roll-off, for ground elevation angles 0deg to 90deg (boresight).

% 0 - 30 degrees should have no gain because no signal

fun = @(x) Gain(x, 120);

ground_el_vals = 0:0.1:90; % degrees
plot(ground_el_vals, fun(ground_el_vals));
title("Roll-Off Loss vs Ground Elevation");
xlabel("Ground Elevation (Degrees)");
ylabel("Gain (dB)");

%% Part B
% Combine your results in part (a) with the impact of changing the range
% between the satellite and the ground station to show the total change in
% received power due to beam off-roll and changing range for elevantion
% angles between 0 and 90

fun = @(x) Gain(x, 120);
fun2 = @(x) SpaceLoss(Range(x,700));
y = fun(ground_el_vals) + fun2(ground_el_vals);
plot(ground_el_vals, y);
title("Roll-Off and Space Loss vs Ground Elevation");
xlabel("Ground Elevation (Degrees)");
ylabel("Gain (dB)");

function gain = Gain(ground_el, full_angle)
    alpha = full_angle;
    theta0 = alpha/2;
    
    theta = 90 - ground_el;
    if theta > theta0
        gain = 0;
    else
        gain = -12.*(theta/alpha).^2;
    end
end

function range = Range(ground_el, h)
    rE = 6378; % km
    range = sqrt((rE + h)^2 - (rE*cosd(ground_el)).^2) - rE*sind(ground_el);
end

function space_loss = SpaceLoss(Range)
    space_loss = 10*log10(1./(4*pi*Range.^2));
end







