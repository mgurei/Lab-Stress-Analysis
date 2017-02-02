% Generate heart rate variability timeine from QRS time locations
% Author: George Tanev
% Input:
%  - qrs_time: timeline of detected QRS locations
%  - f_resample: resampling frequency to have a uniform time between HRV
%  samples
%  - fs: sampling frequency of ECG signal
% Output:
%  - HRV: computed from difference of time locations of QRS complexes 
%  (interbeat intervals) 
%  - HRV: resampled HRV interpolated by f_resample to have uniform time 
%  between HRV samples
%  - duration: duration of HRV timeline
function [HRV qrs_loc HRV_resample qrs_loc_resample] = s_get_HRV(qrs_time, f_resample, fs)

    % Timeline of ECG
    t=1/fs:1/fs:length(qrs_time)*1/fs;

    % Find time instances of QRS complex
    qrs_loc=t(qrs_time~=0);
    HRV=diff(qrs_loc); %Calculate interbeat differences
    
    % TASK: Resample HRV
    % Reampling HRV for frequency analysis
    qrs_loc_resample=0;  % Resampled time series
    HRV_resample=0;     % Resamples HRV
    
end