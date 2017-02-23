function [ newqrs ] = s_detect_filtering(segment, qrs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Get average


size = length(segment); % get the size of the segment for looping
meanv = 0;
n = 0;

for i = 1:size
    if qrs(i) == 1
        meanv = meanv + segment(i);  %if a qrs is detected, add the value of the peak to a value
        n = n + 1;
    end
end
   
meanv = meanv/n;  %Calculate the mean of all detected peaks

lim = meanv*0.85;  %Calculate a limit for the peaks to be discarded

for i = 1:size  % Loop through the segment again
    if qrs(i) == 1
        if segment(i) < lim  % if a qrs is below the limit, remove it
            qrs(i) = 0;
        end 
    end
end
newqrs = qrs;
end

