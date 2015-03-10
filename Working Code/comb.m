function [linfilt Mwave boxfilt] = comb(in,M)
%Function Description
%Takes in an input matrix of EMG (in) and the number of previous frames to
%consider (M) and returns the filtered response as a single vector
%(linfilt), the Mwave, and the filtered response as an array (boxfilt). in
%needs to be in array format with rows being successive frames of data
%between stimulation pulses and columns being successive timepoints within
%a frame.

%This isn't technically a comb filter. It's really just an average of
%previous frames.


%Parameters
[m n]=size(in);

%Initializing
linfilt=[];
boxfilt=[];
Mwave=[];

%Function Begin
for i=1:m-M
    
    temp=mean(in(i:i+M-1,:));
    
    linfilt(end+1:end+n)=in(i+M,:)-temp;
    
    boxfilt(i,:)=in(i+M,:)-temp; %this was previously just temp.
    
    Mwave(end+1:end+n)=temp;
    
end