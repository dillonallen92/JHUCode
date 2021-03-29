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
sigmaVec_1 = zeros(10);
sigmaVec_2 = zeros(10);


 sigmaVec_1 = AllanVariance(clock1.tClock, clock1.tVec, clock1.dt, tauVec);
 sigmaVec_2 = AllanVariance(clock2.tClock, clock2.tVec, clock2.dt, tauVec);


%% Part 3
% Plot the Allan Deviation using the data computed above for each clock. 
% Use the figure in the homework file to determine what type of clock
% these data files represent

% plot(log10(tauVec), sigmaVec_1);
plot(log10(tauVec), sigmaVec_2);


%% Functions
% Allan Variance Function
function sigma_y = AllanVariance(clockVec, timeVec,dt, tauVec)
    % clockVec - Clock time from the structure given
    % timeVec - true time vector given
    % dt - timestep for given matrix
    driftVec = clockVec - timeVec;
    for j = 1 : length(tauVec)
        tau = tauVec(j);
        n = round(tau);
        count = length(driftVec)/n;
        sum = 0;
        for i = 1 : count-2
            sum = sum + (driftVec(i+2*n) - 2*driftVec(i+n) + driftVec(i))^2;
        end
        sigmaSquared = (1/(2*(count-2)*tau^2))*sum;
        sigma_y(j) = log10(sqrt(sigmaSquared));
    end
end

