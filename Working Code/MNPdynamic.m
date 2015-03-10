% function [framesOfTime framesOfVolEMG framesOfEMG framesOfmnpEMG istimTime] = MNPdynamic(tincr,U,Nstim,plotlevel,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
function [t mnpEMG framesOfEMG framesOfmnpEMG stimClock volemg mwaves] = MNPdynamic(tincr,U,Nstim,plotlevel,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% model of force, emg produced by dynamic excitation

freqstim = 25;
ISI = .04;
iISI = round(ISI/tincr);

tstart = 0.1; % start of simulted contraction, s
tstop = 3; % end of simulation, s
tstim = 1.5;
istim = floor(tstim/tincr) + 1;
ivol=[floor((tstim-0.6)/tincr) floor((tstim-0.1)/tincr)];
icomb=[floor((tstop-0.5)/tincr) floor((tstop)/tincr)];
t = 0:tincr:tstop;
totsamp = length(t);
coefficientOfVariation = .05;
stimT = t(istim:end) - t(istim);

istimTime = 1:iISI:totsamp;
stimClock = zeros(1,totsamp);
stimClock(istimTime) = 1;       % This clock is used to synchronize the stimulus pulses              

% stimTime = 0:ISI:tstop;
% istimTime = ceil(stimTime/tincr) + 1;
% stimClock = zeros(1,totsamp);
% stimClock(istimTime) = 1;

%% shape of voluntary excitation
switch U
    case 'exp'
        excitation = exp(t*log(max(exc)+1)/tstop) - 1;
    case 'ramp'
        excitation = t*max(exc)/tstop;
    otherwise
        excitation = zeros(1,totsamp) + U*max(exc);
end
%% shape of stimulation
switch Nstim
    case 'exp'
        stim = floor(exp(stimT*log(N+1)/stimT(end)) - 1);
    case 'ramp'
        stim = floor(stimT*N/stimT(end));
    otherwise
        stim = zeros(1,length(stimT)) + Nstim;
end

stim = [zeros(1,totsamp-length(stim)) stim];
%% generate force vs. time and EMG vs. time vectors for each active motor unit
smuforce = zeros([N totsamp]);
smuemg = zeros([N totsamp]);
smuMwaves = zeros([N totsamp]);

[~,leng2] = size(smuap);
gnorm=(1-exp(-2*(0.4)^3))/0.4;      % normalized force gain
stimTimesArray = zeros(N,1);

for k = 1:N         % enter loop for each motor neuron
    texc = [0];     % vector of action potential times
    stimTimes = [0];    % vector of action potential times caused by stimulation
    volTimes = [0];   
    twitch=exp(1-t/Tc(k));  % twitch force
    for j = 1:totsamp-1     % enter loop of each stimulus time
        stimulated(1:N)=false; % initialize all motor neurons to false (not excited by stimulation at this time point)
        if stim(j)>0;
            for i=1:1:stim(j)
                stimulated(stimorder(i))=true;
            end
        end
        deltaexc = excitation(j) - RET(k);      % voluntary excitation - recruitment threshold for this motor neuron
        
        if stimulated(k)        
            if deltaexc > 0 && min(Rmin + m*deltaexc,Rmax(k)) > freqstim        % Case 1: fv > fs
                useStim = 0;
                murate = min(Rmin + m*deltaexc,Rmax(k));
                refract = 1./murate;
                tFire = t(j) + (1/murate)*coefficientOfVariation*randn(1,1);    % adds variation to firing time
                if t(j) - texc(end) >= refract
                    texc = [texc tFire];
                    volTimes = [volTimes tFire];
                end
            elseif stimClock(j) == 1        % Case 2: fs >= fv and stimulation is firing at this time point
                murate = freqstim;
                refract = 1./murate;
                tFire = t(j);
                if t(j) - volTimes(end) >= refract
                    texc = [texc tFire];
                    stimTimes = [stimTimes tFire];
                end
            end
%             keyboard
%             if ~wait && t(j) - texc(end) >= refract
%                 texc = [texc tFire];
%                 if useStim
%                     stimTimes = [stimTimes tFire];
%                 else
%                     volTimes = [volTimes tFire];
%                 end
%             end
        elseif deltaexc > 0             % Case 3: fv = 0
            murate = min(Rmin + m*deltaexc,Rmax(k)); 
            refract = 1./murate;
            tFire = t(j)+ (1/murate)*coefficientOfVariation*randn(1,1);
            
            if t(j) - texc(end) >= refract
                texc = [texc tFire];
                volTimes = [volTimes tFire];
            end
        end
    end
    
    texc = texc(2:end);
    stimTimes = stimTimes(2:end);
    lengthOfRow = length(stimTimes);
    stimTimesArray(k,1:lengthOfRow) = stimTimes;
    volTimes = volTimes(2:end);
    
    for j=1:1:length(texc); %j indexes the stimulus times
        istart=round(texc(j)/tincr);
        
        % determine force gain
        if j==1;
            g=1;
        else
            Tnorm=Tc(k)/(texc(j)-texc(j-1));
            if Tnorm<0.4;
                g=1;
            else
                g=((1-exp(-2*(Tnorm)^3))/Tnorm)/gnorm;
            end
        end
        scale=g*P(k)/Tc(k);
        for i=istart:totsamp; %i indexes the time since the start of the impulse                 
            smuforce(k,i)=smuforce(k,i)+scale*(t(i)-t(istart))*twitch(i-istart+1);
            if (i-istart+1 <= leng2)
                smuemg(k,i) = smuemg(k,i) + smuap(k,i-istart+1);
            end
        end
    end
    
    for j=1:1:length(stimTimes); %j indexes the stimulus times
        istart=round(stimTimes(j)/tincr);
        for i=istart:totsamp; %i indexes the time since the start of the impulse                 
            if (i-istart+1 <= leng2)
                smuMwaves(k,i) = smuMwaves(k,i) + smuap(k,i-istart+1);
            end
        end
    end
end

mnpForce = sum(smuforce,1);
mnpEMG = sum(smuemg,1);
mwaves = sum(smuMwaves,1);

numberOfFrames = floor(length(mnpEMG)/iISI);
% EMGstim = mnpEMG(icomb(1):icomb(2));
% numberOfFramesStim = floor(length(EMGstim)/iISI);

extra = mod(length(mnpEMG),iISI);
framesOfEMG = reshape(mnpEMG(1:end-extra),iISI,numberOfFrames);
framesOfEMG = framesOfEMG';

% extra = mod(length(EMGstim),iISI);
% framesOfEMGStim = reshape(EMGstim(1:end-extra),iISI,numberOfFramesStim);
% framesOfEMGStim = framesOfEMGStim';

if plotlevel >= 1; %plot individual responses
    figure
    subplot(4,1,1)
    plot(t,mnpForce); ylabel('force');
    hold;
    subplot(4,1,2)
    plot(t,mnpEMG); ylabel('EMG');
    subplot(4,1,3)
    plot(t,excitation); ylabel('excitation');
    subplot(4,1,4)
    plot(t,stim); ylabel('stimulation');
    hold;
end

framesOfmnpEMG = reshape(mnpEMG(1:end-extra),iISI,numberOfFrames);
framesOfmnpEMG = framesOfmnpEMG';

volemg = mnpEMG - mwaves;
framesOfVolEMG = reshape(volemg(1:end-extra),iISI,numberOfFrames);
framesOfVolEMG = framesOfVolEMG';

framesOfTime = reshape(t(1:end-extra),iISI,numberOfFrames);
framesOfTime = framesOfTime';
rectvolemg = abs(volemg);

if plotlevel==2;
    figure;
    subplot(3,1,1),plot(t,mwaves);
    ylabel('mwaves');
    subplot(3,1,2),plot(t,volemg);
    ylabel('volemg');
    subplot(3,1,3),plot(t,rectvolemg);
    ylabel('rectvolemg');
    xlabel('time, s');
end


%calculate average values in the two time periods
Fvol=mean(mnpForce(ivol(1):ivol(2)));
Fcomb=mean(mnpForce(icomb(1):icomb(2)));
%use absolute values for EMGs
% Evol = mean(abs(mnpEMG(ivol(1):ivol(2))));
Evolcomb = mean(abs(mnpEMG(icomb(1):icomb(2))));
Evolstim = mean(rectvolemg(icomb(1):icomb(2)));
