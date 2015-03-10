function [outfilt sumbs trash] = sennel(matrix,M)
%Function Description
%Takes in an input matrix (matrix) of emg where the different rows are the different
%emg frames and then removes the Mwave based on the number of frames
%described (M). Outputs are outfilt which is the linear vector output,
%sumbs is the Mwave, and trash isn't anything, but it was supposed to be
%the matrix form of the filtered output.


%Parameters
[m n]=size(matrix); %m is frames and n is time points

%Initializing
bigphi=[];
phirs=[];
sumbs=zeros(1,n);

%Function 

matrix=flipud(matrix);

for r=1:M+1 %changed from 0:M
    
    for s=2:M+1 %changed from 1:M %since the first row is out actual frame of interest and there is never a 
        
        temp=[];
        
        for i=1:n
            
            temp(end+1)=matrix(r,i)*matrix(s,i);
            
        end
        
        bigphi(r,s-1) = sum(temp);%s-1 to make it line up and then we later remove the zero line for bigphi
        
    end
    
end

theta=bigphi(1,:);

bigphi(1,:)=[];

b=bigphi\theta';

outfilt = [];
for j=1:M
    
    sumbs=sumbs+b(j)*matrix(j+1,:);
    outfilt = [outfilt matrix(j,:)-sumbs];
        
end
    
% outfilt=matrix(1,:)-sumbs;

trash=[];
    