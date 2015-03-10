%ModelSetup.m
%Creates the motor neuron pool model and calculates i/o properties
%Saves a file ***.mat to be used by later analysis programs
%
% plotlevel variable controls the number of plots generated
% 0 = no plots; 1 = high-level summary plots; 2 = adds detailed time function plots
%% set up motor neuron pool recruitment and rate modulation model
clear;
plotlevel=2;
model=2; % see MNPinit for model parameters
tincr=0.0005; %time increment, s, in force and EMG summation array
[exc steps murate Rmin Rmax m P Tc RET N stimorder smuap] = MNPinit4(model,tincr, plotlevel);

%% caclulate voluntary static i/o properties
plotlevel=1; %
Exc=[0 0.05:0.05:1.0]; %Exc is fraction of full excitation
Nexc=length(Exc);

MNPdynamic(tincr,Exc(end),0,2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
MNPforce4(tincr,Exc(end),0,2,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);

for i=1:1:Nexc;
    [Fvolonly(i) temp Evolonly(i) temp1]=MNPforce4(tincr, Exc(i),0,plotlevel,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);%nothing stimulated
end
vFmax=Fvolonly(Nexc);%maximal voluntary force
vEmax=Evolonly(Nexc); %maximal voluntary EMG
Fvolnorm=Fvolonly/vFmax; %force normalized to vFmax
Evolnorm=Evolonly/vEmax; %EMG normalized to vEmax
Uvalues=[0 0.05:0.05:1];%desired equally spaced force increments
Excnorm=interp1(Fvolnorm, Exc,Uvalues,'spline'); % this finds the excitations
                                        %required for 5% increments of vFmax
figure;
hold;
plot(Exc,Fvolonly,'ro');
xlabel('excitation');ylabel('voluntary force');
figure;
plot(Uvalues,Excnorm);
xlabel('excitation');ylabel('inputs for equal force increments');
%% calculate stimulated static i/o properties
stim=0:5:N;
for i=1:length(stim);
    [temp Fstimonly(i)]=MNPdynamic(tincr, 0,stim(i),2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);%no voluntary excitation
    [temp Fstimonly(i)]=MNPforce4(tincr, 0,stim(i),2,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);%no voluntary excitation
end
sFmax=Fstimonly(length(stim));%save for comparison with vFmax
figure;
hold;
plot(stim,Fstimonly,'ro');
xlabel('number of stimulated axons');ylabel('stimulated force');

%estimate numbers of axons to get increments of 10%vFmax up to 60%vFmax
stimincr=0.1*vFmax;
stimvalues= 0:stimincr:6*stimincr ;
stimnum=floor(interp1(Fstimonly,stim,stimvalues,'spline')); 

save ('model.mat');
% %% single occlusion test for separate vol and stim components
plotlevel=2; %
% try at single levels of excitation and stimulation
[Fvoltemp Fcombtemp Evoltemp Ecombtemp]=MNPforce4(tincr, Excnorm(5),stimnum(5),plotlevel,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);
Fstimtemp=Fcombtemp-Fvoltemp;

% stop
%% occlusion tests
% plotlevel=2; %
% Fvol=zeros(7,7);
% Fcomb=zeros(7,7);
% Fstim=zeros(7,7);
% Evol=zeros(7,7);
% Ecomb=zeros(7,7);
% for j=2:1:7;%0 to 60% of stimulated forces
%     for i=1:1:7;%0 to 60% of voluntary force levels
%         [Fvol(i,j) Fcomb(i,j) Evol(i,j) Ecomb(i,j)]=MNPdynamic(tincr, Excnorm(2*i-1),stimnum(j),plotlevel,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
%         [Fvol(i,j) Fcomb(i,j) Evol(i,j) Ecomb(i,j)]=MNPforce4(tincr, Excnorm(2*i-1),stimnum(j),plotlevel,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);
%         Fstim(i,j)=Fcomb(i,j)-Fvol(i,j);
%     end
% end
% %%
% figure;
% hold;
% symcolor=['yo';'ko';'ro';'go';'bo';'co';'mo'];
% for j=2:1:7
%     subplot(3,2,j-1), plot(Fvol(2:7,j)/vFmax,Fstim(2:7,j)/Fstim(2,j),symcolor(j-1));
%     axis([0 1 0 1.1]);
%     xlabel('norm vol force');ylabel('norm stim incr');
% end
% 
% figure;
% hold;
% symcolor=['ko-';'ro-';'yo-';'go-';'bo-';'co-';'mo-'];
% for j=2:1:7
%     vFnorm(:,j)=Fvol(2:7,j)/vFmax;
%     sFnorm(:,j)=Fstim(2:7,j)/Fstim(1,j);
%     plot(vFnorm(:,j),sFnorm(:,j),symcolor(j-1,:));
% end
% axis([0 0.75 0 1.1]);
% xlabel('normalized voluntary force');ylabel('normalized stimulated force increment');
% [s errmsg]=sprintf('Model %d',model); %
% title(s);
% %% trial plots
% figure %force vs. emg
% plot(Evol/vEmax,Fvol/vFmax);
% hold;
% plot(Evolnorm,Fvolnorm);
% axis([0 1.1 0 1.1]);
% xlabel('normalized voluntary EMG');ylabel('normalized voluntary force');
% 
