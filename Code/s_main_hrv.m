%-----------------------Loading ECG Data----------------------------------%
% Version: Student
% This file loads the ECG data and the starting points of every
% stage for each subject (see paper on Acute Stress by Tanev et all)
close all; clear all; clc;
%% Load Data
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
classes={'neutral';'negative';'mental task';'neutral';'baseline'};
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
ecg=ecg_data(1,~isnan(ecg_data(1,:)));

Segment{1} = ecg(1:startPoint(1,1));
Segment{2} = ecg(startPoint(1,1)+1:startPoint(1,2));
Segment{3} = ecg(startPoint(1,2)+1:startPoint(1,3));
Segment{4} = ecg(startPoint(1,3)+1:startPoint(1,4));
Segment{5} = ecg(startPoint(1,4)+1:length(ecg));

% Old Algorithm
% for i = 1:5
%     qrs{i} = s_detect_qrs(Segment{i},b_low,b_high,b_avg,delay);
%     qrs{i} = s_detect_filtering(Segment{i},qrs{i});
% end

% New Algorithm
for i=1:5
    [qrs_amp_raw,qrsi{i},delay]=pan_tompkin(Segment{i},fs,0);
    numNN{i} = length(qrsi{i});
    qrs{i} = transforMat(qrsi{i},length(Segment{i}));
end

figure;
plot(1:length(qrs{1}),5000*qrs{1},1:length(Segment{1}),Segment{1});

% line([startPoint(1,1),startPoint(1,1)],[0,5000],'LineWidth',1,'Color','Green');
% line([startPoint(1,2),startPoint(1,2)],[0,5000],'LineWidth',1,'Color','Green');
% line([startPoint(1,3),startPoint(1,3)],[0,5000],'LineWidth',1,'Color','Green');
% line([startPoint(1,4),startPoint(1,4)],[0,5000],'LineWidth',1,'Color','Green');
% 
% title('ECG with detected QRS')
% 
% % TASK: SEGMENT ECG APPROPRIATELY AND FIND QRS FOR EACH INDIVIDUAL SEGMENT
% % ECG SHOULD BE SEGMENTED ACCORDING TO THE EXPERIMENTAL STAGES
% 
%% Part 2: Calculating HRV
% Example of computing HRV
f_resample=8;

% TASK: FOR EACH SEGMENT OF THE QRS DETECTED SEGMENTS CALCULATE THE HRV 
% TACHOGRAM AND RESAMPLE IT EVENLY TO 8 Hz 
% (make changes in s_get_HRV use the interp1 function to interpolate the HRV)

% This version of the s_get_HRV online returns correct values for HRV and
% qrs_loc_time. You need to make appropriate changes to reutrn HRV_resample and
% qrs_loc_time_resample.

% HRV and HRV resampled calculation
for i = 1:5
[~, ~, HRV_resample{i}, qrs_loc_time_resample{i}] = s_get_HRV(qrs{i}, f_resample, fs);
    HRV_resample{i}(isnan(HRV_resample{i})) = [];
end
% % Plotting HRV
% figure;
% plot(qrs_loc_time(2:end), HRV);
% title('HRV tachogram');
% xlabel('Time of RR (seconds)')
% xlabel('RR interval (seconds)');
% 
% % Plotting HRV resampled
% figure;
% plot(qrs_loc_time_resample, HRV_resample);
% title('HRV resampled tachogram');
% xlabel('Time of RR (seconds)')
% xlabel('RR interval (seconds)');

% THIS FOURIER TRANSFORM WILL ONLY WORK WHEN THE APPROPRIATE CHANGES WERE
% MADE TO THE s_get_HRV function FOR THE HRV_resample and
% qrs_loc_time_resample

figure; hold on
for i = 1:5
temp=HRV_resample{i}(~isnan(HRV_resample{i}));
HRV_psd{i}=abs((fftshift(fft(temp)).^2)/length(temp)); %PSD (s^2/Hz) estimation
freq_axis=-f_resample/2:f_resample/(length(HRV_psd{i})):...
        (f_resample/2-f_resample/(length(HRV_psd{i})));
plot(freq_axis,HRV_psd{i});
end


legend(classes);    
title('Power Spectra distribution')
xlabel('Frequency (Hz)'); ylabel('PSD (s^2/Hz)');
axis([0,0.4,0,20]);
hold off

%% Part 3: Extracting HRV Parameters
% 
% TASK: EXTRACT PARAMETERS FOR EACH OF THE HRV SEGMENTS (FOR EACH
% EXPERIMENTAL STAGE). TIME DOMAIN, FREQUENCY DOMAIN AND NON-LINEAR 
% PARAMETERS (OPTIONAL) 

% % Time Domain Parameters
% % mean
% % standard deviation of NN intervals
% % RMS of difference between adjacent NN intervals
% % STD of difference between adjacent NN intervals
% % NN50 - number of NN intervals that differ by less than 50 ms
% % pNN50 - percentage of NN intervals that differ by less than 50 ms from
% % all NN intervals
limit_50ms = 0.05;

for i = 1:5
HRV_mean{i} = mean(HRV_resample{i});
HRV_std{i} = std(HRV_resample{i});

    for j = 1:(length(HRV_resample{i})-1)
       NN_diff{i}(j) = HRV_resample{i}(j+1) - HRV_resample{i}(j);      
    end
NN_diff_rms{i} = rms(NN_diff{i});
NN_diff_std{i} = std(NN_diff{i});
NN50{i} = HRV_resample{i} (HRV_resample{i} < limit_50ms);
pNN50{i} = 100* NN50{i}/ numNN{i};
end
% % Frequency Domain Parameters
% % Total HRV Positive frequency range
% % VLF - HRV from 0.00 to 0.04
% % LF - HRV from 0.05 to 0.15 Hz (normalize)
% % HF - HRV from 0.16 to 0.40 Hz (normalize)
% Ratio of LF to HF
tot_range = [0 f_resample/2-0.01];
VLF_range = [0 0.04];
LF_range = [0.05 0.15];
HF_range = [0.16 0.4];

for i = 1:5
tot_band_power{i} = bandpower(HRV_psd{i},f_resample, tot_range);
VLF_band_power{i} = bandpower(HRV_psd{i},f_resample, VLF_range);
LF_band_power{i} = bandpower(HRV_psd{i},f_resample, LF_range);
HF_band_power{i} = bandpower(HRV_psd{i},f_resample, HF_range);
ratioLH{i} = LF_band_power{i}/HF_band_power{i};
end



% %% Part 4: Present HRV results

% %% NON-LINEAR PARAMETERS
for i = 1:5  % Approximate entropy
    AppEn{i} = ApEn( 2, 0.2*SDNN{i}, Segment{i}); 
    %   dim : embedded dimension
    %   r : tolerance (typically 0.2 * std)
    %   data : time-series data
    % The smaller the value, the more regular and preditive the sequence is
end

for i = 1:5  % Sample entropy
    SaEn{i} = SampEn( 2, 0.2*SDNN{i}, Segment{i}); 
    %   dim     : embedded dimension
    %   r       : tolerance (typically 0.2 * std)
    %   data    : time-series data
end

for i = 1:5  % Correlation dimention
    [D2, Cm, epsilon] = corrdim(Segment{i},2,tau,epsilon,sloperange); 
    %   m - embedding dimension
    %   tau - delay time lag
    
end



[Alpha1 Alpha2]=DFA_main(DATA); % Detrended fluctuation analysis
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
%   A is the alpha in the paper
%   D is the dimension of the time series
%   n can be changed to your interest




%% Part 4: Present HRV results

% 
% % TASK: PRESENT HRV PARAMETERS IN BOX AND WHISKER PLOTS BETWEEN THE FOUR
% % EXPERIMENTAL STAGES. YOU HAVE THE OPTION TO ATTEMPT BAYESION
% % CLASSIFICATION OF THE RESULTS.
