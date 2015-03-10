function [force,EMG,mWaves] = runSim(ISI,...
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
                                     maxFiringRates)

% Check that all necessary variables are present
if ismember(0,ismember({'ISI',...
                        'timeIncrement',...
                        'duration',...
                        'volitionalSignal',...
                        'stimulatedSignal',...
                        'stimulationStart',...
                        'recruitmentThresholds',...
                        'slopeOfFiringRate',...
                        'numberOfMotorNeurons',...
                        'sMUAP',...
                        'contractionTimes',...
                        'thresholdFiringRate',...
                        'coefficientOfVariation',...
                        'peakTwitchForces',...
                        'maxFiringRates'},who))
    error('Insufficient parameters available')
end

iISI = round(ISI/timeIncrement);
time = 0:timeIncrement:duration;
totalSamples = length(time);
istimulationStart = floor(stimulationStart/timeIncrement) + 1;
stimulationTime = time(istimulationStart:end) - time(istimulationStart);

stimulusTimes = 1:iISI:totalSamples;
stimulationClock = zeros(1,totalSamples);
stimulationClock(stimulusTimes) = 1;

%% shape of voluntary excitation
switch volitionalSignal
    case 'exp'
        vSignal = exp(time*log(101)/duration) - 1;
    case 'ramp'
        vSignal = time*100/duration;
    otherwise
        vSignal = zeros(1,totalSamples) + volitionalSignal*100;
end

%% shape of stimulated excitation
switch stimulatedSignal
    case 'exp'
        sSignal = floor(exp(stimulationTime*log(numberOfMotorNeurons + 1)/stimulationTime(end)) - 1);
    case 'ramp'
        sSignal = floor(stimulationTime*numberOfMotorNeurons/stimulationTime(end));
    otherwise
        sSignal = zeros(1,length(stimulationTime)) + stimulatedSignal;
end

sSignal = [zeros(1,totalSamples - length(sSignal)) sSignal];

%% generate firing time vectors

motorUnitForce = zeros([numberOfMotorNeurons totalSamples]);
motorUnitEMG = zeros([numberOfMotorNeurons totalSamples]);
motorUnitMwaves = zeros([numberOfMotorNeurons totalSamples]);

stimulationOrder = randperm(numberOfMotorNeurons);
normalizedForceGain = (1 - exp(-2*0.4^3))/0.4;
[~,sMUAPlength] = size(sMUAP);

for sample = 1:totalSamples - 1
    isStimulated = zeros(totalSamples,numberOfMotorNeurons);
    if sSignal(sample) > 0
        for index = sSignal(sample)
            isStimulated(sample,stimulationOrder(index)) = true;
        end
    end
end
    
for motorNeuron = 1:numberOfMotorNeurons
    excitationTimes = 0;
    stimulationTimes = 0;
    
    twitchForce = exp(1 - time/contractionTimes(motorNeuron));
    
    for sample = 1:totalSamples - 1
        
        degreeOfExcitation = vSignal(sample) - recruitmentThresholds(motorNeuron);
        possibleFiringRate = min(thresholdFiringRate + slopeOfFiringRate*degreeOfExcitation,maxFiringRates(motorNeuron));
        
        if isStimulated(sample,motorNeuron)
            if degreeOfExcitation > 0 && possibleFiringRate > 1/ISI
                firingRate = possibleFiringRate;
                refractoryPeriod = 1./firingRate;
                firingTime = time(sample) + refractoryPeriod*coefficientOfVariation*randn(1,1);
                if time(sample) - excitationTimes(end) >= refractoryPeriod
                    excitationTimes = [excitationTimes firingTime];
                end
            elseif stimulationClock(sample) == 1
                refractoryPeriod = ISI;
                firingTime = time(sample);
                if time(sample) - excitationTimes(end) >= refractoryPeriod
                    excitationTimes = [excitationTimes firingTime];
                    stimulationTimes = [stimulationTimes firingTime];
                end
            end
        elseif degreeOfExcitation > 0 
            firingRate = possibleFiringRate;
            refractoryPeriod = 1./firingRate;
            firingTime = time(sample) + refractoryPeriod*coefficientOfVariation*randn(1,1);
            if time(sample) - excitationTimes(end) >= refractoryPeriod
                    excitationTimes = [excitationTimes firingTime];
            end
        end
    end
    
    excitationTimes = excitationTimes(2:end);
    stimulationTimes = stimulationTimes(2:end);
    
    for pulseTime = 1:length(excitationTimes) %j indexes the stimulus times
        startingIndex = round(excitationTimes(pulseTime)/timeIncrement);
        if pulseTime == 1
            forceGain = 1;
        else
            normalizedContractionTime = contractionTimes(motorNeuron)/(excitationTimes(pulseTime) - excitationTimes(pulseTime - 1));
            if normalizedContractionTime < 0.4;
                forceGain = 1;
            else
                forceGain = ((1 - exp(-2*(normalizedContractionTime)^3))/normalizedContractionTime)/normalizedForceGain;
            end
        end
        scalingFactor = forceGain*peakTwitchForces(motorNeuron)/contractionTimes(motorNeuron);
        for index = startingIndex:totalSamples %i indexes the time since the start of the impulse                 
            motorUnitForce(motorNeuron,index) = motorUnitForce(motorNeuron,index) + scalingFactor*(time(index) - time(startingIndex))*twitchForce(index - startingIndex + 1);
            if (index - startingIndex + 1 <= sMUAPlength)
                motorUnitEMG(motorNeuron,index) = motorUnitEMG(motorNeuron,index) + sMUAP(motorNeuron,index - startingIndex + 1);
            end
        end
    end
    
%     for pulseTime = 1:length(excitationTimes)
%         
%         startingIndex = round(excitationTimes(pulseTime)/timeIncrement);
%         
%         if pulseTime == 1
%             forceGain = 1;
%         else
%             normalizedContractionTime = contractionTimes(motorNeuron)/(excitationTimes(pulseTime) - excitationTimes(sample - 1));
%             if normalizedContractionTime < 0.4
%                 forceGain = 1;
%             else
%                 forceGain = ((1 - exp(-2*(normalizedContractionTime)^3))/normalizedContractionTime)/normalizedForceGain;
%             end
%         end
%         
%         scalingFactor = forceGain*peakTwitchForces(motorNeuron)/contractionTimes(motorNeuron);
%         
%         for index = startingIndex:totalSamples
%             
%             motorUnitForce(motorNeuron,index) = motorUnitForce(motorNeuron,index) + scalingFactor*(time(index) - time(startingIndex))*twitchForce(index - startingIndex + 1);
%             
%             if (index - startingIndex + 1 <= sMUAPlength)
%                 motorUnitEMG(motorNeuron,index) = motorUnitEMG(motorNeuron,index) + sMUAP(motorNeuron,index - startingIndex + 1);
%             end
%         end
%     end
    
    for pulseTime = 1:length(stimulationTimes)
        
        startingIndex = round(stimulationTimes(pulseTime)/timeIncrement);
        
        for index = startingIndex:totalSamples
            if (index - startingIndex + 1 <= sMUAPlength)
                motorUnitMwaves(motorNeuron,index) = motorUnitMwaves(motorNeuron,index) + sMUAP(motorNeuron,index - startingIndex + 1);
            end
        end
    end
end

force = sum(motorUnitForce,1);
EMG = sum(motorUnitEMG,1);
mWaves = sum(motorUnitMwaves,1);