%% model parameters
numberOfMotorNeurons = 120; % number of motor units
rangeOfThresholds = 40; % range of recruitment thresholds (stated nominal value = 40% in paper)
rangeOfTwitchAmplitude = 100; % range of peak twitch amplitudes
rangeOfTwitchContractionTime = 3; % range of twitch contraction times
longestContractionTime = 0.09; % longest motor unit contraction time, s
shortestContractionTime = longestContractionTime/rangeOfThresholds; % shortest contraction time
maxFRofHighest = 52; % max firing rate of highest threshold unit
thresholdFiringRate = 8; % all units have the same threshold firing rate

%% simulation parameters
timeIncrement = .0004; % sampling interval
ISI = .04;
duration = 3;
stimulationStart = 1.5;

volitionalSignal = .4;
stimulatedSignal = 40;

%% MN pool

a = log(rangeOfThresholds)/numberOfMotorNeurons;
b = log(rangeOfTwitchAmplitude)/numberOfMotorNeurons;
c = log(rangeOfTwitchAmplitude)/log(rangeOfTwitchContractionTime);

% all units have the same slope of firing rate versus excitation
slopeOfFiringRate = (maxFRofHighest - thresholdFiringRate)/(100 - rangeOfThresholds);

i = 1:numberOfMotorNeurons;
recruitmentThresholds = exp(a*i);
peakTwitchForces = exp(b*i);
contractionTimes = longestContractionTime*(1./peakTwitchForces).^(1/c);
maxFiringRates = 1.5./contractionTimes;

% single motor unit action potentials
MUAPduration = 0.012;
lambda = MUAPduration/5;
timeSteps = -MUAPduration/2:timeIncrement:MUAPduration/2;
[~,APlength] = size(timeSteps);
sMUAP = zeros(numberOfMotorNeurons,APlength);
normalizedMUAP = timeSteps.*exp(-(timeSteps.^2)/lambda^2);

for motorNeuron = 1:numberOfMotorNeurons
    sMUAP(motorNeuron,:) = peakTwitchForces(motorNeuron)*normalizedMUAP(:);
end

coefficientOfVariation = .05;

[force, EMG, mWaves] = runSim(ISI,...
      timeIncrement,...
      duration,...
      volitionalSignal,...
      stimulatedSignal,...
      stimulationStart,...
      recruitmentThresholds,...
      slopeOfFiringRate,...
      numberOfMotorNeurons,...
      sMUAP,...
      contractionTimes,...
      thresholdFiringRate,...
      coefficientOfVariation,...
      peakTwitchForces,...
      maxFiringRates);