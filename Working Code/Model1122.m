%%
clear;
plotlevel=0;
model=2; % see MNPinit for model parameters
tincr=.0004; %time increment, s, in force and EMG summation array
[exc steps murate Rmin Rmax m P Tc RET N stimorder smuap] = MNPinit4(model,tincr, plotlevel);

%%
[time,EMG,~,~,clock,volemg,mwaves] = MNPdynamic(tincr,.1,'ramp',2,exc,steps,Rmin,Rmax,m,P,Tc,RET,N,smuap,stimorder);
% figure; plot(time,volemg,'r'); xlim([1.3 end]); ylim([-1 1]); title('Volitional EMG'); xlabel('Time (s)'); ylabel('Potential')
% figure; plot(time,mwaves,'g'); xlim([1.3 end]); ylim([-1 1]); title('M-waves'); xlabel('Time (s)'); ylabel('Potential')
% figure; plot(time,EMG,'k'); xlim([1.5 1.7]); ylim([-.5 .5]); title('Combined EMG'); xlabel('Time (s)'); ylabel('Potential')

[frames extra] = adjustFrames(EMG,clock);


% 
% toFix = [58:61 64:66 69:73 76 78:83 89:94 100:105];
% for i = toFix
%     temp = frames(i,:);
%     frames(i,:) = [temp(2:end) 0];
% end

% keyboard

M = 3;

[filteredVector Mwave filteredArray] = SimpleGSFilt(frames,M);

timeShiftedFilteredArray = padarray(filteredArray,[M 0],NaN,'pre');

[timeFrames] = adjustFrames(time,clock);
[volFrames] = adjustFrames(volemg,clock);


% adjustedFiltered = [timeShiftedFilteredArray extra];
% 
[rows columns] = size(timeShiftedFilteredArray);

EMGvector = reshape(timeShiftedFilteredArray',1,rows*columns);
timeVector = reshape(timeFrames',1,rows*columns);
EMGadjusted = reshape(frames',1,rows*columns);
volemg = reshape(volFrames',1,rows*columns);

figure
plot(timeVector,EMGadjusted,'LineWidth',3)
hold on
plot(timeVector,EMGvector,'g','LineWidth',2)
plot(timeVector,volemg,'r','LineWidth',1)
hold on

residual = EMGvector - volemg;
rms(M) = sqrt(nanmean(residual.^2))
