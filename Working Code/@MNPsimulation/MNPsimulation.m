classdef MNPsimulation
    % MNPsimulation properties:
    %   duration - duration of the simulation, in seconds
    %   samplingRate - time increment between samples, in seconds
    %   stimulationFrequency - frequency of stimulation pulses, in Hertz
    %   voluntaryWaveform - value or shape of the volitional excitation, in
    %   percent of maximum excitation
    %   stimulationWaveform - value or shape of the stimulation , in
    %   percent of number of motor neurons
    %   numberOfMotorNeurons - integer number of motor units in the motor
    %   unit pool
    %   rangeOfThresholds - range of recruitment thresholds
    %   rangeOfTwitchAmplitude - range of peak twitch amplitudes
    %   rangeOfTwitchContractionTime - range of twitch contraction times
    %   longestContractionTime - longest motor unit contraction time, in
    %   seconds
    %   maxFRofHighest - max firing rate of highest threshold unit, in
    %   Hertz
    %   thresholdFiringRate - initial firing rate of all motor units at
    %   threshold, in Hertz
    %   plotlevel - parameter to determine depth of plots to display,
    %   integer from 1 to 3
    properties
        duration = 3;
        stimulationStartTime = 1.5;
        samplingRate = 2000;
        stimulationFrequency = 35;
        voluntaryLevel = .1;
        stimulationLevel = .5;
        numberOfMotorNeurons = 120;
        rangeOfThresholds = 40;
        rangeOfTwitchAmplitude = 100;
        rangeOfTwitchContractionTime = 3;
        longestContractionTime = .09;
        maxFRofHighest = 52;
        thresholdFiringRate = 8;
        maxExcitation = 100;
        coefficientOfVariance = .05;
        MUAPduration = .012;
        plotLevel = 2;
    end
    
    % dependent properties
    properties (SetAccess = 'private')
        timeIncrement
        time
        numberOfSamples
        slopeOfFiringRate
        stimulationWaveform
        voluntaryWaveform
        recruitmentThresholds
        peakTwitchForces
        contractionTimes
        maxFiringRates
        stimulationPulses
        stimulationOrder
        MUAPs
    end
    
    % methods to keep user-controlled parameters within bounds
    methods (Static)
        function isSinglePositiveDouble = isSinglePositiveDouble(number)
            if isequal(size(number),[1 1]) && number > 0 && isequal(class(number),'double')
                isSinglePositiveDouble = 1;
            else
                error('Property value must be a positive double of size [1 1]')
            end
        end
        
        function isFraction = isFraction(number)
            if MNPsimulation.isSinglePostiveDouble(number) && number <= 1
                isFraction = 1;
            else
                error('Property value must be between 0 and 1')
            end
        end
        
        function isInteger = isInteger(number)
            if MNPsimulation.isSinglePositiveDouble(number) && mod(number,1) == 0
                isInteger = 1;
            else
                error('Property value must be an integer')
            end
        end
    end
    
    % set methods for user-controlled parameters
    methods
        function obj = set.duration(obj,duration)
            if MNPsimulation.isSinglePositiveDouble(duration)
                obj.duration = duration;
            end
        end
        
        function obj = set.samplingRate(obj,samplingRate)
            if MNPsimulation.isSinglePositiveDouble(samplingRate)
                obj.samplingRate = samplingRate;
            end
        end
        
        function obj = set.voluntaryLevel(obj,voluntaryWaveform)
            if MNPsimulation.isFraction(voluntaryWaveform)
                obj.voluntaryLevel = voluntaryWaveform;
            end
        end
        
        function obj = set.stimulationLevel(obj,stimulationWaveform)
            if MNPsimulation.isFraction(stimulationWaveform)
                obj.stimulationLevel = stimulationWaveform;
            end
        end
        
        function obj = set.numberOfMotorNeurons(obj,numberOfMotorNeurons)
            if MNPsimulation.isInteger(numberOfMotorNeurons)
                obj.numberOfMotorNeurons = numberOfMotorNeurons;
            end
        end
        
        function obj = set.rangeOfThresholds(obj,rangeOfThresholds)
            if MNPsimulation.isSinglePositiveDouble(rangeOfThresholds)
                obj.rangeOfThresholds = rangeOfThresholds;
            end
        end
        
        function obj = set.rangeOfTwitchAmplitude(obj,rangeOfTwitchAmplitude)
            if MNPsimulation.isSinglePositiveDouble(rangeOfTwitchAmplitude)
                obj.rangeOfTwitchAmplitude = rangeOfTwitchAmplitude;
            end
        end
        
        function obj = set.rangeOfTwitchContractionTime(obj,rangeOfTwitchContractionTime)
            if MNPsimulation.isSinglePositiveDouble(rangeOfTwitchContractionTime)
                obj.rangeOfTwitchContractionTime = rangeOfTwitchContractionTime;
            end
        end
        
        function obj = set.longestContractionTime(obj,longestContractionTime)
            if MNPsimulation.isSinglePositiveDouble(longestContractionTime)
                obj.longestContractionTime = longestContractionTime;
            end
        end
        
        function obj = set.maxFRofHighest(obj,maxFRofHighest)
            if MNPsimulation.isSinglePositiveDouble(maxFRofHighest)
                obj.maxFRofHighest = maxFRofHighest;
            end
        end
        
        function obj = set.thresholdFiringRate(obj,thresholdFiringRate)
            if MNPsimulation.isSinglePositiveDouble(thresholdFiringRate)
                obj.thresholdFiringRate = thresholdFiringRate;
            end
        end
        
        function obj = set.maxExcitation(obj,maxExcitation)
            if MNPsimulation.isSinglePositiveDouble(maxExcitation)
                obj.maxExcitation = maxExcitation;
            end
        end
        
        function obj = set.coefficientOfVariance(obj,coefficientOfVariance)
            if MNPsimulation.isFraction(coefficientOfVariance)
                obj.coefficientOfVariance = coefficientOfVariance;
            end
        end
        
        function obj = set.stimulationStartTime(obj,stimulationStartTime)
            if MNPsimulation.isSinglePositiveDouble(stimulationStartTime)
                obj.stimulationStartTime = stimulationStartTime;
            end
        end
        
        function obj = set.MUAPduration(obj,MUAPduration)
            if MNPsimulation.isSinglePositiveDouble(MUAPduration)
                obj.MUAPduration = MUAPduration;
            end
        end
        
        function obj = set.plotLevel(obj,plotLevel)
            if MNPsimulation.isInteger(plotLevel) && plotLevel < 4
                obj.plotLevel = plotLevel;
            end
        end
    end
    
    % defining dependent properties
    methods
        function timeIncrement = get.timeIncrement(obj)
            timeIncrement = 1/obj.samplingRate;
        end
        
        function time = get.time(obj)
            time = 0:obj.timeIncrement:obj.duration;
        end
        
        function slopeOfFiringRate = get.slopeOfFiringRate(obj)
            slopeOfFiringRate = (obj.maxFRofHighest - obj.thresholdFiringRate)/(obj.maxExcitation - obj.rangeOfThresholds);
        end
        
        function numberOfSamples = get.numberOfSamples(obj)
            numberOfSamples = length(obj.time);
        end
        
        function stimulationWaveform = get.stimulationWaveform(obj)
            if MNPsimulation.isSinglePositiveDouble(obj.stimulationLevel)
                stimulationWaveform = ones(1,obj.numberOfSamples)*obj.stimulationLevel;
            else
                stimulationWaveform = obj.stimulationLevel;
            end
        end

        function voluntaryWaveform = get.voluntaryWaveform(obj)
            if MNPsimulation.isSinglePositiveDouble(obj.voluntaryLevel)
                voluntaryWaveform = ones(1,obj.numberOfSamples)*obj.voluntaryLevel;
            else
                voluntaryWaveform = obj.voluntaryLevel;
            end
        end
        
        function numberStimulated = numberStimulated(obj,sample)
            numberStimulated = floor(obj.stimulationWaveform(sample)*obj.numberOfMotorNeurons);
        end
        
        function maxFiringRates = get.maxFiringRates(obj)
            maxFiringRates = 1.5./obj.contractionTimes;
        end
        
        function recruitmentThresholds = get.recruitmentThresholds(obj)
            a = log(obj.rangeOfThresholds)/obj.numberOfMotorNeurons;
            recruitmentThresholds = exp(a*(1:obj.numberOfMotorNeurons));
        end
        
        function peakTwitchForces = get.peakTwitchForces(obj)
            peakTwitchForces = 1:obj.numberOfMotorNeurons;
            b = log(obj.rangeOfTwitchAmplitude)/obj.numberOfMotorNeurons;
            peakTwitchForces = exp(b*peakTwitchForces);
        end
        
        function contractionTimes = get.contractionTimes(obj)
            c = log(obj.rangeOfTwitchAmplitude)/log(obj.rangeOfTwitchContractionTime);
            contractionTimes = obj.longestContractionTime*(1./obj.peakTwitchForces).^(1/c);
        end
        
        function stimulationPulses = get.stimulationPulses(obj)
            stimulationTime = 0:1/obj.stimulationFrequency:obj.duration;
            stimulationTimeIndecies = ceil(stimulationTime/obj.timeIncrement) + 1;
            stimulationPulses = zeros(1,obj.numberOfSamples);
            stimulationPulses(stimulationTimeIndecies) = 1;
        end
        
        function stimulationOrder = get.stimulationOrder(obj)
            stimulationOrder = randperm(obj.numberOfMotorNeurons);
        end
        
        function MUAPs = get.MUAPs(obj)
            lambda = obj.MUAPduration/5;
            timeSteps = -obj.MUAPduration/2:obj.timeIncrement:obj.MUAPduration/2;
            [~, MUAPlength] = size(timeSteps);
            MUAPs = zeros(obj.numberOfMotorNeurons,MUAPlength);
            normalizedMUAP = timeSteps.*exp(-(timeSteps.^2)/lambda^2);
            
            for motorNeuron = 1:obj.numberOfMotorNeurons
                MUAPs(motorNeuron,:) = obj.peakTwitchForces(motorNeuron)*normalizedMUAP;
            end
        end
    end
    
    methods
        function [EMGvector MwaveVector forceVector] = evaluate(obj)
            time = obj.time;
            timeIncrement = obj.timeIncrement;
            stimulationWaveform = obj.stimulationWaveform;
            voluntaryWaveform = obj.voluntaryWaveform;
            slopeOfFiringRate = obj.slopeOfFiringRate;
            numberOfSamples = obj.numberOfSamples;
            recruitmentThresholds = obj.recruitmentThresholds;
            contractionTimes = obj.contractionTimes;
            peakTwitchForces = obj.peakTwitchForces;
            normalizedForceGain = (1-exp(-2*(0.4)^3))/0.4;
            MUAPs = obj.MUAPs;
            stimulationOrder = obj.stimulationOrder;
            singleUnitForce = zeros(obj.numberOfMotorNeurons,numberOfSamples);
            singleUnitEMG = zeros(obj.numberOfMotorNeurons,numberOfSamples);
            singleUnitMwaves = zeros(obj.numberOfMotorNeurons,numberOfSamples);
            [~,MUAPlength] = size(MUAPs);
            
            for motorNeuron = 1:obj.numberOfMotorNeurons
                excitationTimes = 0;
                stimulationTimes = 0;
                voluntaryTimes = 0;
                
                twitchForce=exp(1 - time/contractionTimes(motorNeuron)); %#ok<*PROP>
                
                for sample = 1:numberOfSamples - 1
                    isStimulated = zeros(1,obj.numberOfMotorNeurons); %initialize them all to false
                    if stimulationWaveform(sample) > 0;
                        for i = 1:obj.numberStimulated(sample)
                            isStimulated(stimulationOrder(i)) = true;
                        end
                    end
                    deltaExcitation = voluntaryWaveform(sample)*obj.maxExcitation - recruitmentThresholds(motorNeuron);

                    if isStimulated(motorNeuron)
                        
                        if deltaExcitation > 0 && min(obj.thresholdFiringRate + slopeOfFiringRate*deltaExcitation,obj.maxFiringRates(motorNeuron)) > obj.stimulationFrequency
                            firingRate = min(obj.thresholdFiringRate + slopeOfFiringRate*deltaExcitation,obj.maxFiringRates(motorNeuron));
                            refractoryPeriod = 1/firingRate;
                            firingTime = time(sample) + (refractoryPeriod)*obj.coefficientOfVariance*randn(1,1);
                            if time(sample) - excitationTimes(end) >= refractoryPeriod
                                excitationTimes = [excitationTimes firingTime];
                                voluntaryTimes = [voluntaryTimes firingTime];
                            end
                        elseif obj.stimulationPulses(sample) == 1
                            firingRate = obj.stimulationFrequency;
                            refractoryPeriod = 1/firingRate;
                            firingTime = time(sample);
                            if time(sample) - voluntaryTimes(end) >= refractoryPeriod
                                excitationTimes = [excitationTimes firingTime];
                                stimulationTimes = [stimulationTimes firingTime];
                            end
                        end
            
                    elseif deltaExcitation > 0
                        firingRate = min(obj.thresholdFiringRate + slopeOfFiringRate*deltaExcitation,obj.maxFiringRates(motorNeuron));
                        refractoryPeriod = 1/firingRate;
                        firingTime = time(sample)+ (1/firingRate)*obj.coefficientOfVariance*randn(1,1);

                        if time(sample) - excitationTimes(end) >= refractoryPeriod
                            excitationTimes = [excitationTimes firingTime];
                            voluntaryTimes = [voluntaryTimes firingTime];
                        end
                    end
                end

                excitationTimes = excitationTimes(2:end);
                stimulationTimes = stimulationTimes(2:end);
                voluntaryTimes = voluntaryTimes(2:end);

                for sample = 1:length(excitationTimes); %j indexes the stimulus times
                    startOfMUAP=floor(excitationTimes(sample)/timeIncrement);
                    if sample == 1;
                        forceGain = 1;
                    else
                        normalizedContractionTime=contractionTimes(motorNeuron)/(excitationTimes(sample) - excitationTimes(sample - 1));
                        if normalizedContractionTime < 0.4;
                            forceGain = 1;
                        else
                            forceGain = ((1 - exp(-2*(normalizedContractionTime)^3))/normalizedContractionTime)/normalizedForceGain;
                        end
                    end
                    scaleForce = forceGain*peakTwitchForces(motorNeuron)/contractionTimes(motorNeuron);
                    for i = startOfMUAP:numberOfSamples; %i indexes the time since the start of the impulse                 
                        singleUnitForce(motorNeuron,i) = singleUnitForce(motorNeuron,i) + scaleForce*(time(i) - time(startOfMUAP))*twitchForce(i - startOfMUAP + 1);
                        if (i - startOfMUAP + 1 <= MUAPlength)
                            singleUnitEMG(motorNeuron,i) = singleUnitEMG(motorNeuron,i) + MUAPs(motorNeuron,i - startOfMUAP + 1);
                        end
                    end
                end

                for sample = 1:length(stimulationTimes); %j indexes the stimulus times
                    startOfMUAP=floor(stimulationTimes(sample)/timeIncrement);
                    for i = startOfMUAP:numberOfSamples; %i indexes the time since the start of the impulse                 
                        if (i - startOfMUAP + 1 <= MUAPlength)
                            singleUnitMwaves(motorNeuron,i) = singleUnitMwaves(motorNeuron,i) + MUAPs(motorNeuron,i - startOfMUAP + 1);
                        end
                    end
                end
            end
            forceVector = sum(singleUnitForce,1);
            EMGvector = sum(singleUnitEMG,1);
            MwaveVector = sum(singleUnitMwaves,1);
            
            
        end
    end
end