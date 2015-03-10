function [outfilt Mwave box] = SimpleGSFilt(in,M)
%Function Description
%GS Filter inputs are in-array with emg that includes the M-wave, M-the
%number of previous frames to observe when removing the M-wave, and
%L-the number of samples in a frame. the output outfilt is the 'raw'
%voluntary signal with the M-wave subtracted. Mwave is the removed m-wave,
%and box is the array version of outfilt

%Initializing Parameters
sepdata=in;

%initializing the matrices for the epsilon matrix and the weights matrix.
w=[];
eps=[];
box=[];

%this is the GS Filter stuff. 
data=[]; %initializes the vector that contains the processed emg data
[q r]=size(sepdata); %this tells us how many rows of data there are. each row is one stimulus period.
fullM=[];
Mwaveset=[];
%this was included to account for the order that signals are added in.
iset=M+1:-1:1;

%Function
for count=0:q-M-1 %this way we go through each row except for those that we can't create full matrices for
    
    for i=1:M+1%for i=0:M this does this the first time
        eps(i,1,:)=sepdata(iset(i),:); %initializing the epsilon matrix with raw data
    end    
    
    for m=1:M+1 %this is the same as 0:M-1
        
        for i=1:M-m+1
            temp1(1,:)=eps(i,m,:);
            temp2(:,1)=eps(M-m+2,m,:);
            
            w(m,i)=temp1*temp2/(mag(eps(M-m+2,m,:))^2);
            w;
            eps(i,m+1,:)=eps(i,m,:)-w(m,i)*eps(M-m+2,m,:);
 %           Mset=w(m,i)*eps(M-m+2,m,:); %not the M wave, I need the original portion-minus the raw
        end
        
    end
    
%    Mwaveset(end+1:end+L)=Mset;%this isn't the M wave.
    yo(1,:)=eps(1,M+1,:);%this may need to be M+1
    yol=length(yo);
    data(end+1:end+yol)=yo;
    tempdiff=sepdata(M+1,:)-yo;
    tdiffl=length(tempdiff);
    Mwaveset(end+1:end+tdiffl)=tempdiff;
    sepdata(1,:)=[];
    box(end+1,:)=yo;
end
outfilt=data;
Mwave=Mwaveset;

end

function ans=mag(X)

ans=sqrt(sum((X.^2)));

end

        
        