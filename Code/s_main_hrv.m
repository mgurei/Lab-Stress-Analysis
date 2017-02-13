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
subject = 1;
ecg=ecg_data(subject,~isnan(ecg_data(subject,:)));

Segment{1} = ecg(1:startPoint(subject,1));
Segment{2} = ecg(startPoint(subject,1)+1:startPoint(subject,2));
Segment{3} = ecg(startPoint(subject,2)+1:startPoint(subject,3));
Segment{4} = ecg(startPoint(subject,3)+1:startPoint(subject,4));
Segment{5} = ecg(startPoint(subject,4)+1:length(ecg));

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

% figure;
% plot(1:length(qrs{1}),5000*qrs{1},1:length(Segment{1}),Segment{1});
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
% Her the Heart Rate Variability (HRV) is calculated from the know QRS
% position. The HRV is resampled at 8 Hz to have a uniform time axis.
f_resample=8;

% HRV and HRV resampled calculation
HRV_resample = cell(1,length(classes));
qrs_loc_time_resample = cell(1,length(classes));
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
HRV_psd = cell(1,length(classes));
freq_axis = cell(1,length(classes));

figure; hold on
for i = 1:5
temp=HRV_resample{i}(~isnan(HRV_resample{i}));
HRV_psd{i}=abs((fftshift(fft(temp)).^2)/length(temp)); %PSD (s^2/Hz) estimation
freq_axis{i}=-f_resample/2:f_resample/(length(HRV_psd{i})):...
        (f_resample/2-f_resample/(length(HRV_psd{i})));
plot(freq_axis{i},HRV_psd{i});
end


legend(classes);    
title('Power Spectra distribution')
xlabel('Frequency (Hz)'); ylabel('PSD (s^2/Hz)');
axis([0,0.4,0,20]);
hold off

%% Part 3: Extracting HRV Parameters
% Here time, frequency and non-linear features of the HRV signal are
% calculated. 

% Time Domain Parameters
% mean
% standard deviation of NN intervals
% RMS of difference between adjacent NN intervals
% STD of difference between adjacent NN intervals
% NN50 - number of NN intervals that differ by less than 50 ms
% pNN50 - percentage of NN intervals that differ by less than 50 ms from
% all NN intervals
limit_50ms = 0.05;
NN_mean = cell(1,length(classes));
NN_std = cell(1,length(classes));
NN_diff_rms = cell(1,length(classes));
NN_diff_std = cell(1,length(classes));
NN50 = cell(1,length(classes));
pNN50 = cell(1,length(classes));

for i = 1:5
NN_diff = diff(HRV_resample{i}); % difference between adjancent NN intervals    

HRV_mean{i} = mean(HRV_resample{i}); %mean of HRV
NN_std{i} = std(HRV_resample{i});   %standard deviantion HRV
NN_diff_rms{i} = rms(NN_diff);       %RMS of NN diff
NN_diff_std{i} = std(NN_diff);       %standard deviation od NN diff
NN50{i} = HRV_resample{i} (HRV_resample{i} < limit_50ms);   %NN less the 50ms
pNN50{i} = 100* NN50{i}/ numNN{i};   % percentage NN50 in all HRV
end

% Frequency Domain Parameters
% Total HRV Positive frequency range
% VLF - HRV from 0.00 to 0.04
% LF - HRV from 0.05 to 0.15 Hz (normalize)
% HF - HRV from 0.16 to 0.40 Hz (normalize)
% Ratio of LF to HF
tot_band_power = cell(1,length(classes));
VLF_band_power = cell(1,length(classes));
LF_band_power = cell(1,length(classes));
HF_band_power = cell(1,length(classes));
ratioLH = cell(1,length(classes));
pVLF = cell(1,length(classes));
pLF = cell(1,length(classes));
pHF = cell(1,length(classes));

for i = 1:5
% frenquency used
tot_range = find(freq_axis{i}>=0.00 & freq_axis{i}<f_resample/2);
VLF_range = find(freq_axis{i}>=0.00 & freq_axis{i}<=0.04);
LF_range = find(freq_axis{i}>=0.05 & freq_axis{i}<=0.15);
HF_range = find(freq_axis{i}>=0.16 & freq_axis{i}<=0.4);

%calculations
tot_band_power{i} = (HRV_psd{i}(tot_range)'*HRV_psd{i}(tot_range) )/ (length(tot_range));
VLF_band_power{i} = (HRV_psd{i}(VLF_range)'*HRV_psd{i}(VLF_range) )/ (length(VLF_range));
LF_band_power{i} = (HRV_psd{i}(LF_range)'*HRV_psd{i}(LF_range) )/ (length(LF_range));
HF_band_power{i} = (HRV_psd{i}(HF_range)'*HRV_psd{i}(HF_range) )/ (length(HF_range));
ratioLH{i} = LF_band_power{i}/HF_band_power{i};
% percentage of the bands compared to the total
pVLF{i} = VLF_band_power{i}/tot_band_power{i} * 100;
pLF{i} = LF_band_power{i}/tot_band_power{i} * 100;
pHF{i} = HF_band_power{i}/tot_band_power{i} * 100;

end

% Non-linear parameters calculated here are:
% -AppEn = approximate entropy
% -SaEn = sample entropy 
% -
% - 
dim = 2; % embedded dimension
tau = 0.1; %delay time lag
AppEn = cell(1,length(classes));
SaEn = cell(1,length(classes));

for i = 1:5  
    toll = .2 * NN_std{i}; %tollerance 
    AppEn{i} = ApEn( dim, toll, HRV_resample{i}); % Approximate entropy
    SaEn{i} = SampEn( dim, toll, HRV_resample{i}); %Sample Entropy
%    [D2, Cm, epsilon] = corrdim(HRV_resample{i},dim,tau,epsilon,sloperange); 
%       m - embedding dimension
%       tau - delay time lag

   [Alpha1{i} Alpha2{i}]=DFA_main(HRV_resample{i}); % Detrended fluctuation analysis
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
%   A is the alpha in the paper
%   D is the dimension of the time series
%   n can be changed to your interest

end


%% Part 4: Present HRV results


% TASK: PRESENT HRV PARAMETERS IN BOX AND WHISKER PLOTS BETWEEN THE FOUR
% EXPERIMENTAL STAGES. YOU HAVE THE OPTION TO ATTEMPT BAYESION
% CLASSIFICATION OF THE RESULTS.
