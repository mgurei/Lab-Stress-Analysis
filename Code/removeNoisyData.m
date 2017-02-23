function [ fixedSegment ] = removeNoisyData( Segment )
%UNTITLED Summary of this function goes here
%   A problem that arosed when using the Pan-Thompkins algorithm was if
%   there were  high power noise in the data, the dynamic threshold of the signal would
%   be pushed so high that the real peaks would get ignored for being under
%   the threshold. This was solved by finding the average of the signals
%   and cutting every peak that was high over the average of the signal. 
    
meanV = mean(Segment); %Find 
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