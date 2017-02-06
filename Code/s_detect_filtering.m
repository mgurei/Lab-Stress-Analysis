function [ newqrs ] = s_detect_filtering(segment, qrs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Get average


size = length(segment);
meanv = 0;
n = 0;

for i = 1:size
    if qrs(i) == 1
        meanv = meanv + segment(i);
        n = n + 1;
    end
end
   
meanv = meanv/n;

lim = meanv*0.85;  

for i = 1:size
    if qrs(i) == 1
        if segment(i) < lim
            qrs(i) = 0;
        end 
    end
end
newqrs = qrs;
end

