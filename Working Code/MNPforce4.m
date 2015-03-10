function [Fvol Fcomb Evol Evolcomb Evolstim] = MNPforce4(tincr,U,Nstim,plotlevel,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);
% model of force produced by motor neuron pool
% modified to do separate calculation of vol and stim force, emg
% modified to include calculation and return of EMG - mnpforce3
% modified to use specified order of stimulation recruitment
% based on model published by Fuglevand Winter & Patla 1993 and Zhou and Rymer 2007
% U is the neural input as a fraction of the maximal excitation
% Nstim is the absolute number of motor units recruited by the electrical stimulaltion
%               assume recruitment in reverse order
% stimulation parameters
active=false(N,1); % initialize all motor units as inactive
refract=0.003; %refractory period of stimulation, s
freqstim=35; % stimulation frequency, Hz
ISI=1/freqstim;
%
cv=0.05; %variability of voluntary motor unit firing
% timing parameters
tstart=0.1;% start of simulated contraction, s
tstop=3; %end of simulation, s
tstim=1.5; %start of superimposed electrical stimulation, s
istim=floor(tstim/tincr)+1; %index to start of stimulation
ivol=[floor((tstim-0.6)/tincr) floor((tstim-0.1)/tincr)];%to caclulate average voluntary force
icomb=[floor((tstop-0.5)/tincr) floor((tstop)/tincr)];%to calculate average combined force
t=0:tincr:tstop; %vector of sample times
totsamp=length(t);
%
%% simulate the voluntary and stimulated force response of the MN pool 
excitation = U*max(exc); %fraction of the MVC exc
mnpforce = zeros(1,length(t));% initialize total force vector to zeros
mnpEMG = zeros(1,length(t)); % initialize total EMG vector to zeros
%each row is the force produced by one motor unit at an instant of time
k=1; %k will index motor units
excindex=1;
% find index to excitation array
while exc(excindex)<excitation;
    excindex=excindex+1;
end
%% generate logical vector of whether axons are stimulated (true) or not
%(false)
stimulated(1:N)=false; %initialize them all to false
if Nstim>0;
    for i=1:1:Nstim;
    stimulated(stimorder(i))=true;
    end
end

%% generate force vs. time and EMG vs. time vectors for each active motor unit
smuforce = zeros([N totsamp]);% initialize single motor unit force summing array to zeros
smuemg = zeros([N,totsamp]); %initialize single motor unit EMG summing array to zeros
%calculate force responses to excitations
[~,leng2] = size(smuap);
gnorm=(1-exp(-2*(0.4)^3))/0.4; %for MU force gain normalization
mufreq=zeros(N); %initialize all mu frequencies to zero
for k=1:N;
    %calculate array of precomputed twitch forces
    twitch=exp(1-t/Tc(k));
    %start by creating vectors of times for motor unit excitation
    if RET(k)<excitation; %this section deals with the voluntary group of MUs
        mufreq(k) = murate(excindex,k);
        active(k) = true;%keeps track of which motor units are active
        texcconst = tstart:1/mufreq(k):tstop; %vector of constant stimulus times
        %add ISI variance, use randn to get normal distribution
        texc=texcconst+(1/mufreq(k))*cv*randn(1,length(texcconst));
        %if unit is electrically stimulated, add the stimulation time
        % points, but only if the stimulus period (IPI) is smaller than the motor unit
        % firing period (see Crago 1973).
        if stimulated(k); % handle the voluntary and stimulated group
            active(k)=true;
%            if ISI<1/mufreq(k); %compare vol and stim IPIs 
            if freqstim>mufreq(k); %only substitute stimulus pulses if stim freq exceeds mufreq
                h=1;
                    while texc(h)<=tstim-refract;%skip to the start of stim period
                                                 % and do not stimulate within
                                                 % refractory period of a
                                                 % natural action potential
                        h=h+1;
                    end
                    if texc(h-1)<tstim;
                        texc=[texc(1:h-1) tstim:ISI:tstop];%fill remainder with stim pulses
                    else
                        texc=[texc(1:h) (tstim+ISI):ISI:tstop];
                    end
            end
        end
    elseif stimulated(k); % end of the voluntary loop, start the stimulation only loop
        active(k)=true;
        texc=tstim:ISI:tstop;
    end
    %calculate individual motor unit forces and emgs
    if active(k) ; %only calculate forces of active motor units
        for j=1:1:length(texc); %j indexes the stimulus times
            istart=floor(texc(j)/tincr);
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
    else 
        %no force or emg for inactive motor units
        smuforce(k,:) = 0;
        smuemg(k,:) = 0;
    end
end
    %calculate total force and emg produced by the whole motor neuron pool
mnpforce = sum(smuforce,1);
mnpEMG = sum(smuemg,1);
if plotlevel==2; %plot individual responses
    figure
    subplot(2,1,1)
    plot(t,mnpforce);
    hold;
    subplot(2,1,2)
    plot(t,mnpEMG);
    hold;
end
% Sum the EMGs of the stimulated units
mwaves=zeros(1,totsamp);
for k=1:N;
    if (stimulated(k) && (freqstim>mufreq(k)));
        mwaves(istim:totsamp)=mwaves(istim:totsamp)+smuemg(k,(istim:totsamp));
    end
end
volemg=mnpEMG-mwaves;
rectvolemg=abs(volemg);
%
if plotlevel==2;
    figure;
    subplot(3,1,1),plot(mwaves);
    ylabel('mwaves');
    subplot(3,1,2),plot(volemg);
    ylabel('volemg');
    subplot(3,1,3),plot(rectvolemg);
    ylabel('rectvolemg');
    xlabel('time, s');
end
%calculate average values in the two time periods
Fvol=mean(mnpforce(ivol(1):ivol(2)));
Fcomb=mean(mnpforce(icomb(1):icomb(2)));
%use absolute values for EMGs
Evol = mean(abs(mnpEMG(ivol(1):ivol(2))));
Evolcomb = mean(abs(mnpEMG(icomb(1):icomb(2))));
Evolstim = mean(rectvolemg(icomb(1):icomb(2)));
