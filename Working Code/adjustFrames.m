function [frames extra] = adjustFrames(vector,stimClock)

row = 0;
column = 1;

for index = 1:length(vector)
    if stimClock(index) == 1
        row = row + 1;
        column = 1;
    end
    frames(row,column) = vector(index);
    column = column + 1;
end
% 
% start = 1;
% 
% for i = 1:length(stimClock)
%     for j = start:stimClock(i)
%         frames(row,column) = vector(j);
%         
%         column = column + 1;
%     end
%     start = stimClock(i) + 1;
%     row = row + 1;
%     column = 1;
%     
% end

% minColumns = min(diff(stimClock));

extra = frames(:,end);
% frames = frames(:,1:end-1);