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
classes={'Normal activity';'Negative';'Mental task';'Neutral';'Baseline'};
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

%%     init variables
allNN_diff_rms = zeros(10,4);
allNN_diff_std = zeros(10,4);
allNN_mean = zeros(10,4);
allNN_std = zeros(10,4);
allNN50 = zeros(10,4);
allpNN50 = zeros(10,4);
alltot_band_power = zeros(10,4);
allVLF_band_power = zeros(10,4);
allLF_band_power = zeros(10,4);
allHF_band_power = zeros(10,4);
allratioLH = zeros(10,4);
allpVLF = zeros(10,4);
allpLF = zeros(10,4);
allpHF = zeros(10,4);
allAppEn = zeros(10,4);
allSaEn = zeros(10,4);
allAlpha1 = zeros(10,4);
allAlpha2 = zeros(10,4);


for j=1:10
    %% TEST
    subject = j;
    ecg=ecg_data(subject,~isnan(ecg_data(subject,:)));
    
    s = 10; %seconds of new segment to ignore
    ignore = s*fs;%numbers of samples in a new segment to ignore 
    
    Segment{1} = ecg(1:startPoint(subject,1));
    Segment{2} = ecg(startPoint(subject,1)+1:startPoint(subject,2));
    Segment{3} = ecg(startPoint(subject,2)+1:startPoint(subject,3));
    Segment{4} = ecg(startPoint(subject,3)+1:startPoint(subject,4));
    Segment{5} = ecg(startPoint(subject,4)+1:length(ecg));

    % Old Algorithm
%     for i = 1:5
%         qrs{i} = s_detect_qrs(Segment{i},b_low,b_high,b_avg,delay);
%         qrs{i} = s_detect_filtering(Segment{i},qrs{i});
%     end

    % New Algorithm
    for i=1:5       
        [qrs_amp_raw,qrsi{i},delay]=pan_tompkin(Segment{i},fs,0);
        numNN{i} = length(qrsi{i});
        qrs{i} = transforMat(qrsi{i},length(Segment{i}));
    end

%     figure;hold on;
%     plot((1:length(qrs{1})),5000*qrs{1},'r');
%     plot(startPoint(j,1)+(1:length(qrs{2}))+1,5000*qrs{2},'r');
%     plot(startPoint(j,2)+(1:length(qrs{3}))+1,5000*qrs{3},'r');
%     plot(startPoint(j,3)+(1:length(qrs{4}))+1,5000*qrs{4},'r');
%     plot(startPoint(j,4)+(1:length(qrs{5}))+1,5000*qrs{5},'r');
%     plot(1:length(ecg),ecg);
%     line([startPoint(j,1),startPoint(j,1)],[0,5000],'LineWidth',1,'Color','Green');
%     line([startPoint(j,2),startPoint(j,2)],[0,5000],'LineWidth',1,'Color','Green');
%     line([startPoint(j,3),startPoint(j,3)],[0,5000],'LineWidth',1,'Color','Green');
%     line([startPoint(j,4),startPoint(j,4)],[0,5000],'LineWidth',1,'Color','Green');
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

%     figure; hold on
    for i = 1:5
    temp=HRV_resample{i}(~isnan(HRV_resample{i}));
    HRV_psd{i}=abs((fftshift(fft(temp)).^2)/length(temp)); %PSD (s^2/Hz) estimation
    freq_axis{i}=-f_resample/2:f_resample/(length(HRV_psd{i})):...
            (f_resample/2-f_resample/(length(HRV_psd{i})));
%     plot(freq_axis{i},HRV_psd{i});
    end


%     legend(classes);    
%     title('Power Spectra distribution')
%     xlabel('Frequency (Hz)'); ylabel('PSD (s^2/Hz)');
%     axis([0,0.4,0,20]);
%     hold off

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
    NN_mean = zeros(1,length(classes));
    NN_std = zeros(1,length(classes));
    NN_diff_rms = zeros(1,length(classes));
    NN_diff_std = zeros(1,length(classes));
    NN50 = zeros(1,length(classes));
    pNN50 = zeros(1,length(classes));

    for i = 1:5
    NN_diff = diff(HRV_resample{i}); % difference between adjancent NN intervals    

    NN_mean(1,i) = mean(HRV_resample{i}); %mean of HRV
    NN_std(1,i) = std(HRV_resample{i});   %standard deviantion HRV
    NN_diff_rms(1,i) = rms(NN_diff);       %RMS of NN diff
    NN_diff_std(1,i) = std(NN_diff);       %standard deviation od NN diff
    %if ~isempty(HRV_resample{i} (HRV_resample{i} < limit_50ms))
    NN50(1,i) = numel(NN_diff (abs(NN_diff) > limit_50ms));   %adjancent NN differening  more than 50ms
    pNN50(1,i) = 100* NN50(1,i)/ length(NN_diff);   % percentage NN50 in all HRV
    %end
    end

    % Frequency Domain Parameters
    % Total HRV Positive frequency range
    % VLF - HRV from 0.00 to 0.04
    % LF - HRV from 0.05 to 0.15 Hz (normalize)
    % HF - HRV from 0.16 to 0.40 Hz (normalize)
    % Ratio of LF to HF
    tot_band_power = zeros(1,length(classes));
    VLF_band_power = zeros(1,length(classes));
    LF_band_power = zeros(1,length(classes));
    HF_band_power = zeros(1,length(classes));
    ratioLH = zeros(1,length(classes));
    pVLF = zeros(1,length(classes));
    pLF = zeros(1,length(classes));
    pHF = zeros(1,length(classes));

    for i = 1:5
    % frenquency used
    tot_range = [0.00, 3.99];
    VLF_range = [0.00, 0.04];
    LF_range = [0.05, 0.15];
    HF_range = [0.16, 0.4];

    %calculations
    tot_band_power(1,i) = bandpower(HRV_psd{i}, f_resample, tot_range);
    VLF_band_power(1,i) = bandpower(HRV_psd{i}, f_resample, VLF_range);
    LF_band_power(1,i) = bandpower(HRV_psd{i}, f_resample, LF_range);
    HF_band_power(1,i) = bandpower(HRV_psd{i}, f_resample, HF_range);
    ratioLH(1,i) = LF_band_power(1,i)/HF_band_power(1,i);
    % percentage of the bands compared to the total
    pVLF(1,i) = VLF_band_power(1,i)/tot_band_power(1,i) * 100;
    pLF(1,i) = LF_band_power(1,i)/tot_band_power(1,i) * 100;
    pHF(1,i) = HF_band_power(1,i)/tot_band_power(1,i) * 100;

    end

    % Non-linear parameters calculated here are:
    % -AppEn = approximate entropy
    % -SaEn = sample entropy 
    % -
    % - 
    dim = 2; % embedded dimension
    tau = 0.1; %delay time lag
    AppEn = zeros(1,length(classes));
    SaEn = zeros(1,length(classes));

    for i = 1:5  
        toll = .2 * NN_std(1,i); %tollerance 
        AppEn(1,i) = ApEn( dim, toll, HRV_resample{i}); % Approximate entropy
        SaEn(1,i) = SampEn( dim, toll, HRV_resample{i}); %Sample Entropy
    %    [D2, Cm, epsilon] = corrdim(HRV_resample{i},dim,tau,epsilon,sloperange); 
    %       m - embedding dimension
    %       tau - delay time lag
       [Alpha1(1,i), Alpha2(1,i)]=DFA_main(HRV_resample{i}); % Detrended fluctuation analysis
    end
    

    
    
    %% Adding to full variables
    allNN_diff_rms(j,:) = NN_diff_rms(2:5);
    allNN_diff_std(j,:) = NN_diff_std(2:5);
    allNN_mean(j,:) = NN_mean(2:5);
    allNN_std(j,:) = NN_std(2:5);
    allNN50(j,:) = NN50(2:5);
    allpNN50(j,:) = pNN50(2:5);
    alltot_band_power(j,:) = tot_band_power(2:5);
    allVLF_band_power(j,:) = VLF_band_power(2:5);
    allLF_band_power(j,:) = LF_band_power(2:5);
    allHF_band_power(j,:) = HF_band_power(2:5);
    allratioLH(j,:) = ratioLH(2:5);
    allpVLF(j,:) = pVLF(2:5);
    allpLF(j,:) = pLF(2:5);
    allpHF(j,:) = pHF(2:5);
    allAppEn(j,:) = AppEn(2:5);
    allSaEn(j,:) = SaEn(2:5);
    allAlpha1(j,:) = Alpha1(2:5);
    allAlpha2(j,:) = Alpha2(2:5);
    
    
    j
    
end

%% Part 4: Present HRV results

% TASK: PRESENT HRV PARAMETERS IN BOX AND WHISKER PLOTS BETWEEN THE FOUR
% EXPERIMENTAL STAGES. YOU HAVE THE OPTION TO ATTEMPT BAYESION
% CLASSIFICATION OF THE RESULTS.

% hrv_prep = [HRV_resample{1}; HRV_resample{2}; HRV_resample{3}; HRV_resample{4}; HRV_resample{5}];
% grp = [zeros(1,length(HRV_resample{1})), ones(1,length(HRV_resample{2})), ...
%     2*ones(1,length(HRV_resample{3})), 3*ones(1,length(HRV_resample{4})), ...
%     4*ones(1,length(HRV_resample{5}))];
% figure()
% boxplot(hrv_prep, grp, 'Labels', classes);
% title('HRV for all segment')

% -- Time domain parameters --
classes = {'Negative';'Mental task';'Neutral';'Baseline'};
figure;
subplot(6,1,1);
boxplot(allNN_diff_rms,'Labels',classes)
title('NN Diff RMS')

subplot(6,1,2);
boxplot(allNN_diff_std,'Labels',classes)
title('NN Diff STD')

subplot(6,1,3);
boxplot(allNN_mean,'Labels',classes)
title('NN Mean')

subplot(6,1,4);
boxplot(allNN_std,'Labels',classes)
title('NN Mean')

subplot(6,1,5);
boxplot(allNN50,'Labels',classes)
title('NN50')

subplot(6,1,6);
boxplot(allpNN50,'Labels',classes)
title('pNN50')

% -- Frequency domain parameters -- 

figure;
subplot(4,1,1)
boxplot(alltot_band_power,'Labels',classes)
title('Total band power')

subplot(4,1,2)
boxplot(allVLF_band_power,'Labels',classes)
title('VLF band power')

subplot(4,1,3)
boxplot(allLF_band_power,'Labels',classes)
title('LF band power')

subplot(4,1,4)
boxplot(allHF_band_power,'Labels',classes)
title('HF band power')

figure;
subplot(4,1,1)
boxplot(allratioLH,'Labels',classes)
title('Ratio LF HF')

subplot(4,1,2)
boxplot(allpVLF,'Labels',classes)
title('Percentage VLF')

subplot(4,1,3)
boxplot(allpLF,'Labels',classes)
title('Percentage LF')

subplot(4,1,4)
boxplot(allpHF,'Labels',classes)
title('Percentage HF')

% -- Nonlinear domain parameters --

figure;
subplot(4,1,1)
boxplot(allAppEn,'Labels',classes)
title('Approximate Entropy')

subplot(4,1,2)
boxplot(allSaEn,'Labels',classes)
title('Sample Entropy')

subplot(4,1,3)
boxplot(allAlpha1,'Labels',classes)
title('Alpha 1')

subplot(4,1,4)
boxplot(allAlpha2,'Labels',classes)
title('Alpah 2')

