%-----------------------Loading ECG Data----------------------------------%
% Version: Student
%-- This file loads the ECG data and the starting points of every
% stage for each subject (see paper on Acute Stress by Tanev et all)

%% Load Data
clear all; close all; clc;
load 's_ecg_10_subjects.mat';

%% Declare relevant data for ECG signals
% Start sample point of the four stages (columns) for every subject (rows)
startPoint=[15576+10*60*512 15576+(10+8)*60*512 15576+(10+8+6)*60*512 15576+(10+8+6+6)*60*512;
        19000+10*60*512 19000+(10+7.75)*60*512 19000+(10+7.75+6)*60*512 19000+(10+7.75+6+6)*60*512;
        8794+10*60*512 8794+(10+7.25)*60*512 8794+(10+7.25+6)*60*512 8794+(10+7.25+6+6)*60*512;
        5065+10*60*512 5065+(10+7.1)*60*512 5065+(10+7.1+6)*60*512 5065+(10+7.1+6+6)*60*512;
        10571+10*60*512 10571+(10+7)*60*512 10571+(10+7+6)*60*512 10571+(10+7+6+6)*60*512;
        7049+10*60*512 7049+(10+7.5)*60*512 7049+(10+7.5+6)*60*512 7049+(10+7.5+6+6)*60*512;
        29681+10*60*512 29681+(10+7.15)*60*512 29681+(10+7.15+6)*60*512 29681+(10+7.15+6+6)*60*512;
        11587+10*60*512 11587+(10+8.5)*60*512 11587+(10+8.5+6)*60*512 11587+(10+8.5+6+6)*60*512;
        8821+10*60*512 8821+(10+7)*60*512 8821+(10+7+6)*60*512 8821+(10+7+6+6)*60*512;
        40000+10*60*512 40000+(10+7.25)*60*512 40000+(10+7.25+6)*60*512 40000+(10+7.25+6+6)*60*512];
% Names of four experimental stages
classes=repmat({'negative';'mental_task';'neutral';'baseline'},10,1);
% Names of files with ECG data from 10 subjects
names={'sub1_29_05_2013'; 'sub2_29_05_2013'; 'sub3_29_05_2013'; 
        'sub4_29_05_2013'; 'sub5_29_05_2013'; 'sub6_29_05_2013';
        'sub7_30_05_2013'; 'sub8_30_05_2013'; 'sub9_31_05_2013';
        'sub10_31_05_2013'};
fs=512; % Sampling frequency
%% Part 1: Detect QRS
% Example of QRS detection 
load s_b_coeff.mat; % See desc_filt variable for filter parameters
load s_ecg_10_subjects.mat;
subject_no = 1;
ecg=ecg_data(1,~isnan(ecg_data(subject_no,:)));

% Division of ecg in 5 segments
segment{1} = ecg(1:startPoint(subject_no,1))';
segment{2} = ecg(startPoint(subject_no,1)+1:startPoint(subject_no,2))';
segment{3} = ecg(startPoint(subject_no,2)+1:startPoint(subject_no,3))';
segment{4} = ecg(startPoint(subject_no,3)+1:startPoint(subject_no,4))';
segment{5} = ecg(startPoint(subject_no,4)+1:end)';

% Detection of QRS segments and plotting

for i = 1:5
    qrs{i} = s_detect_qrs(segment{i}, b_low, b_high, b_avg, delay);
    figure(i);
    plot (1:length(qrs{i}),5000*qrs{i}, 1:length(segment{i}), segment{i});
end




% qrs=s_detect_qrs(ecg,b_low,b_high,b_avg,delay);
% figure;
% plot(1:length(qrs),5000*qrs,1:length(ecg),ecg);
% title('ECG with detected QRS')

% TASK: SEGMENT ECG APPROPRIATELY AND FIND QRS FOR EACH INDIVIDUAL SEGMENT
% ECG SHOULD BE SEGMENTED ACCORDING TO THE EXPERIMENTAL STAGES

%% Part 2: Calculating HRV
% Example of computing HRV
f_resample=8;

% TASK: FOR EACH SEGMENT OF THE QRS DETECTED SEGMENTS CALCULATE THE HRV 
% TACHOGRAM AND RESAMPLE IT EVENLY TO 8 Hz 
% (make changes in s_get_HRV use the interp1 function to interpolate the HRV)

% This version of the s_get_HRV online returns correct values for HRV and
% qrs_loc_time. You need to make appropriate changes to reutrn HRV_resample and
% qrs_loc_time_resample.
[HRV qrs_loc_time HRV_resample qrs_loc_time_resample]=s_get_HRV(qrs, f_resample, fs);

figure;
plot(qrs_loc_time(2:end),HRV);
title('HRV tachogram');
xlabel('Time of RR (seconds)')
xlabel('RR interval (seconds)');

% THIS FOURIER TRANSFORM WILL ONLY WORK WHEN THE APPROPRIATE CHANGES WERE
% MADE TO THE s_get_HRV function FOR THE HRV_resample and
% qrs_loc_time_resample

% temp=HRV_resample(~isnan(HRV_resample));
% HRV_psd=abs((fftshift(fft(temp)).^2)/length(temp)); %PSD (s^2/Hz) estimation
% freq_axis=-f_resample/2:f_resample/(length(HRV_psd)):...
%         (f_resample/2-f_resample/(length(HRV_psd)));
% figure;
% plot(freq_axis,HRV_psd);
% xlabel('Frequency (Hz)'); ylabel('PSD (s^2/Hz)');
% axis([0,0.4,0,20]);

%% Part 3: Extracting HRV Parameters

% TASK: EXTRACT PARAMETERS FOR EACH OF THE HRV SEGMENTS (FOR EACH
% EXPERIMENTAL STAGE). TIME DOMAIN, FREQUENCY DOMAIN AND NON-LINEAR 
% PARAMETERS (OPTIONAL) 

% Time Domain Parameters
% mean
% standard deviation of NN intervals
% RMS of difference between adjacent NN intervals
% STD of difference between adjacent NN intervals
% NN50 - number of NN intervals that differ by less than 50 ms
% pNN50 - percentage of NN intervals that differ by less than 50 ms from
% all NN intervals

% Frequency Domain Parameters
% Total HRV Positive frequency range
% VLF - HRV from 0.00 to 0.04
% LF - HRV from 0.05 to 0.15 Hz (normalize)
% HF - HRV from 0.16 to 0.40 Hz (normalize)
% Ratio of LF to HF

%% Part 4: Present HRV results

% TASK: PRESENT HRV PARAMETERS IN BOX AND WHISKER PLOTS BETWEEN THE FOUR
% EXPERIMENTAL STAGES. YOU HAVE THE OPTION TO ATTEMPT BAYESION
% CLASSIFICATION OF THE RESULTS.
