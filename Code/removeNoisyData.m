function [ fixedSegment ] = removeNoisyData( Segment )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
meanV = mean(Segment);
percentage = 0.5;
upperT = meanV + (meanV*percentage);
lowerT = meanV - (meanV*percentage);

for i = 1:length(Segment)
    if Segment(i) > upperT
        Segment(i) = upperT;
    elseif Segment(i) < lowerT
        Segment(i) = lowerT;
    end
end

fixedSegment = Segment;


end

