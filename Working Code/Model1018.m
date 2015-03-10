%
%Creates the motor neuron pool model and calculates i/o properties
%Saves a file ***.mat to be used by later analysis programs
%
% plotlevel variable controls the number of plots generated
% 0 = no plots; 1 = high-level summary plots; 2 = adds detailed time function plots
%% set up motor neuron pool recruitment and rate modulation model
clear;
plotlevel=0;
model=2; % see MNPinit for model parameters
tincr=0.0004; %time increment, s, in force and EMG summation array
[exc steps murate Rmin Rmax m P Tc RET N stimorder smuap] = MNPinit4(model,tincr, plotlevel);

%% test dynamic inputs

% MNPdynamic(tincr,.5,0,1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% MNPdynamic(tincr,'ramp',0,1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% MNPdynamic(tincr,'exp',0,1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% 
% MNPdynamic(tincr,0,N/2,1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% MNPdynamic(tincr,0,'ramp',1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% MNPdynamic(tincr,0,'exp',1,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% 
% MNPdynamic(tincr,.5,N/2,2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% MNPdynamic(tincr,'ramp',N/2,2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
[framesOfTime framesOfVolEMG framesOfEMG framesOfmnpEMG] = MNPdynamic(tincr,.1,N/2,2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
[rows columns] = size(framesOfVolEMG');
volEMGvector = reshape(framesOfVolEMG',1,rows*columns);

[rows columns] = size(framesOfEMG');
EMGvector = reshape(framesOfEMG',1,rows*columns);

[rows columns] = size(framesOfmnpEMG');
MNPvector = reshape(framesOfmnpEMG',1,rows*columns);

[rows columns] = size(framesOfTime');
t = reshape(framesOfTime',1,rows*columns);

M = 3;

for i = 1:columns
    subtractedEMG(i,:) = framesOfmnpEMG(i,:) - framesOfmnpEMG(53,:);
end

% subtractedEMG = framesOfmnpEMG - framesOfmnpEMG(end,:);
figure; plot(subtractedEMG(53:end,:)')
hold on; plot(framesOfmnpEMG(53,:))
[rows columns] = size(framesOfmnpEMG');
MNPvector = reshape(framesOfmnpEMG',1,rows*columns);

[filteredVector Mwave filteredArray] = SimpleGSFilt(framesOfEMG,M);

timeShiftedFilteredVector = padarray(filteredVector,[0 M*rows],NaN,'pre');
timeShiftedFilteredArray = padarray(filteredArray,[M 0],NaN,'pre');

figure
plot(t,volEMGvector)
hold on
plot(t,timeShiftedFilteredVector,'r')
plot(t,MNPvector,'g')

residual = volEMGvector - timeShiftedFilteredVector;
figure
plot(t,residual)

residualArray = framesOfVolEMG - timeShiftedFilteredArray;
rmsVector = sqrt(nanmean(residualArray.^2,2));

meanTimeVector = mean(framesOfTime,2);

figure; plot(meanTimeVector,rmsVector,'.')

figure; plot(t,volEMGvector,'LineWidth',4)
hold on
plot(t,MNPvector,'r','LineWidth',2)
plot(t,timeShiftedFilteredVector,'g','LineWidth',1)
plot(meanTimeVector,rmsVector,'.k')
legend('voluntary','combined','filtered')

n = 85;

figure; plot(t(n*rows:end),volEMGvector(n*rows:end),'LineWidth',4)
hold on
plot(t(n*rows:end),MNPvector(n*rows:end),'r','LineWidth',2)
plot(t(n*rows:end),timeShiftedFilteredVector(n*rows:end),'g','LineWidth',1)
plot(meanTimeVector(n*rows:end),rmsVector(n*rows:end),'.k')
legend('voluntary','combined','filtered')

%*************************************************************************

[framesOfTime framesOfVolEMG framesOfEMG framesOfmnpEMG stimClock] = MNPdynamic(tincr,.1,'ramp',2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
[rows columns] = size(framesOfVolEMG');
volEMGvector = reshape(framesOfVolEMG',1,rows*columns);

[rows columns] = size(framesOfEMG');
EMGvector = reshape(framesOfEMG',1,rows*columns);

[rows columns] = size(framesOfmnpEMG');
MNPvector = reshape(framesOfmnpEMG',1,rows*columns);

[rows columns] = size(framesOfTime');
t = reshape(framesOfTime',1,rows*columns);

M = 3;

for i = 1:columns
    subtractedEMG(i,:) = framesOfmnpEMG(i,:) - framesOfmnpEMG(53,:);
end

% subtractedEMG = framesOfmnpEMG - framesOfmnpEMG(end,:);
% figure; plot(subtractedEMG(53:end,:)')
hold on; plot(framesOfmnpEMG(53,:))
[rows columns] = size(framesOfmnpEMG');
MNPvector = reshape(framesOfmnpEMG',1,rows*columns);

[filteredVector Mwave filteredArray] = SimpleGSFilt(framesOfEMG,M);

timeShiftedFilteredVector = padarray(filteredVector,[0 M*rows],NaN,'pre');
timeShiftedFilteredArray = padarray(filteredArray,[M 0],NaN,'pre');

figure; plot(t,volEMGvector,'LineWidth',4)
hold on
plot(t,MNPvector,'r','LineWidth',2)
plot(t,timeShiftedFilteredVector,'g','LineWidth',1)
plot(meanTimeVector,rmsVector,'.k')
legend('voluntary','combined','filtered')

n = 100;

figure; plot(t(n*rows:end),volEMGvector(n*rows:end),'LineWidth',4)
hold on
plot(t(n*rows:end),MNPvector(n*rows:end),'r','LineWidth',2)
plot(t(n*rows:end),timeShiftedFilteredVector(n*rows:end),'g','LineWidth',1)
plot(meanTimeVector(n*rows:end),rmsVector(n*rows:end),'.k')
legend('voluntary','combined','filtered')
