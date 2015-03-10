function [exc steps murate Rmin Rmax m P Tc RET N stimorder smuap] = MNPinit3(model, tincr, plotlevel);
% creates models of motor neuron pool
% function sets up all of the arrays needed to calculate i/o properties
% also calculates array of scaled motor unit action potentials for later
%    emg synthesis
% based on models published by Fuglevand Winter & Patla 1993 and Zhou and Rymer 2007
% U is the neural input as a fraction of the maximal excitation
%
% MN pool parameters
if model==2; %Fuglevand Winter & Patla 1993 
    N=120; %Number of motor units
    RR=40; %range of recruitment thresholds (stated nominal value = 40% in paper)
    RP=100; %range of peak twitch amplitudes
    RT=3; %range of twitch contraction times
    Tlong=0.09; %longest motor unit contraction time, s
    Tshort=Tlong/RT; %shortest contraction time
    RmaxHigh=52; %max firing rate of highest threshold unit
    RmaxLow=23; %max firing rate of lowest threshold unit
    Rmin=8; %all units have the same threshold firing rate
    Rdecr=10;% decrement in firing rate for model 1 of Fuglevand
end
if model==3; %Zhou and Rymer model smaller number of motor units
    N=30; %Number of motor units
    RR=40; %range of recruitment thresholds (stated nominal value = 40% in paper)
    RP=100; %range of peak twitch amplitudes
    RT=3; %range of twitch contraction times
    Tlong=0.09; %longest motor unit contraction time, s
    Tshort=Tlong/RT; %shortest contraction time
    RmaxHigh=52; %max firing rate of highest threshold unit
    RmaxLow=23; %max firing rate of lowest threshold unit
    Rmin=8; %all units have the same threshold firing rate
    Rdecr=10;% decrement in firing rate for model 1 of Fuglevand
end
if model==4; %Zhou and Rymer model peak firing rate 13-20 Hz
    N=120; %Number of motor units
    RR=40; %range of recruitment thresholds (stated nominal value = 40% in paper)
    RP=100; %range of peak twitch amplitudes
    RT=3; %range of twitch contraction times
    Tlong=0.09; %longest motor unit contraction time, s
    Tshort=Tlong/RT; %shortest contraction time
    RmaxHigh=20; %max firing rate of highest threshold unit
    RmaxLow=13; %max firing rate of lowest threshold unit
    Rmin=8; %all units have the same threshold firing rate
    Rdecr=10;% decrement in firing rate for model 1 of Fuglevand
end
if model==5; % Zhou and Rymer model - smaller range of twitch amplitudes
    N=120; %Number of motor units
    RR=40; %range of recruitment thresholds (stated nominal value = 40% in paper)
    RP=40; %range of peak twitch amplitudes
    RT=3; %range of twitch contraction times
    Tlong=0.09; %longest motor unit contraction time, s
    Tshort=Tlong/RT; %shortest contraction time
    RmaxHigh=52; %max firing rate of highest threshold unit
    RmaxLow=23; %max firing rate of lowest threshold unit
    Rmin=8; %all units have the same threshold firing rate
    Rdecr=10;% decrement in firing rate for model 1 of Fuglevand
end
if model==6; % Zhou and Rymer model - smaller range of twitch contraction times
    N=120; %Number of motor units
    RR=40; %range of recruitment thresholds (stated nominal value = 40% in paper)
    RP=100; %range of peak twitch amplitudes
    RT=2; %range of twitch contraction times
    Tlong=0.09; %longest motor unit contraction time, s
    Tshort=Tlong/RT; %shortest contraction time
    RmaxHigh=52; %max firing rate of highest threshold unit
    RmaxLow=23; %max firing rate of lowest threshold unit
    Rmin=8; %all units have the same threshold firing rate
    Rdecr=10;% decrement in firing rate for model 1 of Fuglevand
end
%% construct MN pool
a=log(RR)/N;
b=log(RP)/N;
c=log(RP)/log(RT);
RET=1:N; %recruitment excitation threshold of each unit
P=1:N; %peak twitch force of each unit, arbitrary units
Tc=1:N; %contraction time of each unit
Rmax=1:N; %maximum firing rate of each unit
%
m=(RmaxHigh-Rmin)/(100-RR); %all units have the same slope of firing rate versus excitation
%
i=1:N;
RET=exp(a*i);
P=exp(b*i);
Tc=Tlong*(1./P).^(1/c); %
% if ratemodel==1; %model 1 of Fuglevand et al.
%     Rmax=RmaxLow-Rdecr*(RET/RET(N));
% end
%if ratemodel==2; %model 2 of Fuglevand et al., peak rate matched to Tc.
   Rmax=(1.5)./Tc;
%end
%if ratemodel=3;
%Rmax=RmaxHigh-(RmaxHigh-RmaxLow)/(Tc(1)-Tc(N))*(Tc-Tc(N)); %Crago's made up method
%Rmax=RmaxLow+(RmaxHigh-RmaxLow)/(Tc(N)-Tc(1))*(Tc-Tc(N)); %Crago's made up method
%
if plotlevel==1;
    figure; %just plots basic values
    subplot(3,2,1), plot(1:N, Rmax);
    axis([0 N 0 RmaxHigh]);
    xlabel('motor unit number');
    ylabel('Max firing rate, Hz');
    [s errmsg]=sprintf('Model %d',model); %
    title(s);
    subplot(3,2,3), plot(1:N,RET);
    axis([0 N 0 100]);
    xlabel('motor unit number');
    ylabel('Recr Excit Threshold');
    subplot(3,2,5), plot(1:N,P);
    axis([0 N 0 125]);
    ylabel('Peak Force');
    xlabel('motor unit number');
    subplot(3,2,2), plot(1:N,Tc);
    axis([0 N 0 0.1]);
    ylabel('Tcont, s');
    xlabel('motor unit number');
    %calculate total twitch force as a function of excitation
    Ptot=1:N;
    Ptot(1)=P(1);
    for i=2:N;
        Ptot(i)=Ptot(i-1)+P(i);
    end
    subplot(3,2,6), plot(RET,Ptot);
    ylabel('total twitch force');xlabel('Recr Excit Threshold');
    subplot(3,2,4),plot(Tc,P);
    ylabel('peak twitch force');xlabel('Tcont');
    axis([Tshort Tlong 0 100]);
end
%
%calculate firing rate of each unit over the full range of excitation
exc=[RET ceil(RET(N)+.5):100]; %Creates a vector of excitations for each motor unit threshold
%and each integer value of excitation up to 100 after full recruitment
steps=length(exc); % number of excitation levels in exc
murate=zeros(steps,N);%preallocate space for rate vs. excitation, MU number
for j=1:steps;
    for i=1:N; %step through all motor units and 
                      %calculate firing rates at each excitation level
        deltaexc=exc(j)-RET(i);
        if deltaexc<0
            murate(j,i)=0;
        else
            murate(j,i)=min(Rmin+m*deltaexc,Rmax(i));
        end
    end
end
if plotlevel==1;
    figure;
    hold;
    for i=1:N;
        plot(exc(:),murate(:,i));
        xlabel('Excitation');ylabel('MU Rate, Hz');
    end
    [s errmsg]=sprintf('Model %d',model); %
    title(s);
end

stimorder=randperm(N); %randomizes order of recruitment by stimulation
%% generate array (smuap) of scaled single motor unit action potentials
% x=linspace(-4,4,15/(1000*tincr)+2);% time duration of an unscaled     
% impulse = exp(-x.^2/2)/sqrt(2*pi);% undifferentiated impulse
% for i=1:N;
% smuap(i,:)=P(i)*-diff(diff(impulse)); %doubly differentiated
% end
% if plotlevel == 2;
%     figure
%     plot(smuap')
% end
%tsteps=-0.03:tincr:0.03; % s, time steps for calculating the normalized muap
muapduration=0.012; % s, duration of muaps
lambda=muapduration/5; % value 5 determined by trial and error to get correct duration
tsteps=-muapduration/2:tincr:muapduration/2; % s, time steps for calculating the normalized muap
[temp aplength]=size(tsteps);
smuap=zeros(N,aplength);
smuapnorm=tsteps.*exp(-(tsteps.^2)/lambda^2);
for i=1:N;
    smuap(i,:)=P(i)*smuapnorm(:);
end
if plotlevel==2;
    figure;
    plot(tsteps,smuap');
end