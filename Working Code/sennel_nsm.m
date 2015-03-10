function [outfilt Mset outarray] = sennel(matrix,M)
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
outfilt=[];
outarray=[];
Mset=[];
sumbs=zeros(1,n);

%Function 

matrix=flipud(matrix);
fullmatrix=matrix;

for k=1:m-M
    matrix=fullmatrix(k:k+M,:);
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
    k
    b=bigphi\theta'
    
    %outfilt = [];
    sumbs=zeros(1,n);
    for j=1:M
        
        sumbs=sumbs+b(j)*matrix(j+1,:);
        %     outfilt = [outfilt matrix(j,:)-sumbs]; This was in here. I moved it
        %     outside the loop.
        
    end
    
    outfilt = [matrix(j,:)-sumbs outfilt]; %This may be backwards.
    Mset = [sumbs Mset];
    outarray = [outarray;matrix(j,:)-sumbs];
    % outfilt=matrix(1,:)-sumbs;
end
    