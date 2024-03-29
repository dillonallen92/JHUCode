clear, clc, close all;

%% Load the 2 clock data files
% Each file contains timestep (dt), true time vector (tVec) and 
% clock time (tClock)

clock1 = load('HW8_Prob1_clock_1.mat');
clock2 = load('HW8_Prob1_clock_2.mat');

%% Part 2
% Compute the Allan Variance of the clocks for 10 tau values spaced
% evenly in log space 

tauVec = logspace(1,4,10);

 [variance_1, sigmaVec_1] = AllanVariance(clock1.tClock, clock1.tVec, clock1.dt, tauVec);
 [variance_2, sigmaVec_2] = AllanVariance(clock2.tClock, clock2.tVec, clock2.dt, tauVec);


%% Part 3
% Plot the Allan Deviation using the data computed above for each clock. 
% Use the figure in the homework file to determine what type of clock
% these data files represent
subplot(211)
plot(log10(tauVec), sigmaVec_1);
xlabel("log(\tau)");
ylabel("log(\sigma_{y})");
subplot(212)
plot(log10(tauVec), sigmaVec_2);
xlabel("log(\tau)");
ylabel("log(\sigma_{y})");
%% Functions
% Allan Variance Function
function [variance, sigma_y] = AllanVariance(clockVec, timeVec,dt, tauVec)
    % clockVec - Clock time from the structure given
    % timeVec - true time vector given
    % dt - timestep for given matrix
    sigma_y = zeros(length(tauVec));
    driftVec = timeVec - clockVec;
    for j = 1 : length(tauVec)
        tau = tauVec(j);
        n = round(tau/dt);
        count = round(length(driftVec)/n);
        sum = 0;
        for i = 1 : count-2
            sum = sum + (driftVec(i+2*n) - 2*driftVec(i+n) + driftVec(i))^2;
        end
        variance = (1/(2*(count-2)*tau^2))*sum;
        sigma_y(j) = log10(sqrt(variance));
    end
end

%% Discussion
% I believe the first plot is close to passive H Maser while the second one
% is quartz? I think I really messed up this script so I am probably off
% base. 
