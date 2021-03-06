%Analysis.m
%Relies on having a .mat file from a previously executed program called
%ModelSetup.m
% suggested process for running from the command window is to
% first clear all variables, then open the .mat file by double clicking on it
% and then execute analysis.m
%
% plotlevel controls the number of plots generated by the function MNPforce
% 0 = no plots; 1 = high-level summary plots; 2 = adds detailed time function plots
%
%% occlusion tests
plotlevel=2; %
Fvol=zeros(6,6);
Fcomb=zeros(6,6);
Fstim=zeros(6,6);
Evol=zeros(6,6);
Ecomb=zeros(6,6);
for j=2:1:6;%0 to 50% of stimulated forces
    for i=1:1:6;%0 to 50% of voluntary force levels
        [Fvol(i,j) Fcomb(i,j) Evol(i,j) Ecomb(i,j) Evolstim(i,j)]=MNPforce4(tincr, Excnorm(2*i-1),stimnum(j),plotlevel,exc,steps,murate,P,Tc,RET,N,smuap,stimorder);
        Fstim(i,j)=Fcomb(i,j)-Fvol(i,j);
    end
end
figure;
hold;
symcolor=['ko-';'ro-';'yo-';'go-';'bo-';'co-'];
for j=2:1:6;
    vFnorm(:,j)=Fvol(2:6,j)/vFmax;
    sFnorm(:,j)=Fstim(2:6,j)/Fstim(1,j);
    plot(vFnorm(:,j),sFnorm(:,j),symcolor(j-1,:));
end
axis([0 0.75 0 1.1]);
xlabel('normalized voluntary force');ylabel('normalized stimulated force increment');
[s errmsg]=sprintf('Model %d',model); %
title(s);
%% trial plots 
figure %force vs. emg
plot(Evol/vEmax,Fvol/vFmax);
hold;
plot(Evolnorm,Fvolnorm);
axis([0 1.1 0 1.1]);
xlabel('normalized voluntary EMG');ylabel('normalized voluntary force');
%% plot voluntary EMGs as a function of effort and stimulation level
figure;
hold;
symcolor=['ko-';'ro-';'yo-';'go-';'bo-';'co-';'mo-'];
for j=2:1:6
    vEnorm(:,j)=Evol(2:6,j)/vEmax;
    vEstimnorm(:,j)=Evolstim(2:6,j)/vEmax;
    plot(vEnorm(:,j),vEstimnorm(:,j),symcolor(j-1,:));
end
axis([0 .6 0 .6]);
line([0 .6],[0 .6],'Color','k'); %unitary gain line
xlabel('normalized voluntary EMG');ylabel('normalized voluntary EMG during stimulation');
[s errmsg]=sprintf('Model %d',model); %
title(s);

