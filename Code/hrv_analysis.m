%% Modified s_main_hrv for classification testing

function data_out = hrv_analysis(ecg, fs, startPoint, classes)
%% Part 1: Detect QRS

%load specific data and segmentation
Segment{1} = ecg(1:startPoint(1,1));
Segment{2} = ecg(startPoint(1,1)+1:startPoint(1,2));
Segment{3} = ecg(startPoint(1,2)+1:startPoint(1,3));
Segment{4} = ecg(startPoint(1,3)+1:startPoint(1,4));
Segment{5} = ecg(startPoint(1,4)+1:length(ecg));

% Pan Tompkin's QRS algorithm detection
for i=1:5
    [~,qrs_num{i},~] = pan_tompkin(Segment{i},fs,0);
    numNN{i} = length(qrs_num{i});
    qrs{i} = transforMat(qrs_num{i},length(Segment{i}));
end

%% Part 2: Calculating HRV
% Here the Heart Rate Variability (HRV) is calculated from the know QRS
% position. The HRV is resampled at 8 Hz in order to have a uniform time axis.

f_resample=8;
HRV_resample = cell(1,length(classes));
qrs_loc_time_resample = cell(1,length(classes));

for i = 1:5
    [~, ~, HRV_resample{i}, qrs_loc_time_resample{i}] = s_get_HRV(qrs{i}, f_resample, fs);
    HRV_resample{i}(isnan(HRV_resample{i})) = [];
end

% Power spectra density estimation of the HRV (resempled at 8Hz)
HRV_psd = cell(1,length(classes));
freq_axis = cell(1,length(classes));

for i = 1:5
    temp=HRV_resample{i}(~isnan(HRV_resample{i}));
    HRV_psd{i}=abs((fftshift(fft(temp)).^2)/length(temp)); %PSD (s^2/Hz) estimation
    freq_axis{i}=-f_resample/2:f_resample/(length(HRV_psd{i})):...
        (f_resample/2-f_resample/(length(HRV_psd{i})));
end

%% Part 3: Extracting HRV Parameters
% Here time, frequency and non-linear features of the HRV signal are
% calculated.

% Time Domain Parameters:
% NN_mean - mean of NN intervals
% NN_std - standard deviation of NN intervals
% NN_diff_rms - RMS of difference between adjacent NN intervals
% NN_diff_std - STD of difference between adjacent NN intervals
% NN50 - number of NN intervals that differ by more than 50 ms
% pNN50 - percentage of NN intervals that differ by more than 50 ms

limit_50ms = 0.05;
NN_mean = zeros(1,length(classes));
NN_std = zeros(1,length(classes));
NN_diff_rms = zeros(1,length(classes));
NN_diff_std = zeros(1,length(classes));
NN50 = zeros(1,length(classes));
pNN50 = zeros(1,length(classes));
tags = {'NN_mean' 'NN_std' 'NN_rms' 'NN_std' 'NN50' 'pNN50'};

for i = 1:5
    NN_diff = diff(HRV_resample{i}); % difference between adjancent NN intervals
    NN_mean(1,i) = mean(HRV_resample{i}); %mean of HRV
    NN_std(1,i) = std(HRV_resample{i});   %standard deviantion HRV
    NN_diff_rms(1,i) = rms(NN_diff);       %RMS of NN diff
    NN_diff_std(1,i) = std(NN_diff);       %standard deviation od NN diff
    NN50(1,i) = numel(NN_diff (abs(NN_diff) > limit_50ms));   %adjancent NN differening  more than 50ms
    pNN50(1,i) = 100* NN50(1,i)/ length(NN_diff);   % percentage NN50 in all HRV
end

data_out.time = table(NN_mean', NN_std', NN_diff_rms', NN_diff_std', NN50', ...
    pNN50', classes);

% Frequency Domain Parameters
% VLF_band_power - Total HRV Positive frequency range
% VLF_band_power - HRV from 0.00 to 0.04
% LF_band_power - HRV from 0.05 to 0.15 Hz (normalize)
% HF_band_power - HRV from 0.16 to 0.40 Hz (normalize)
% ratioLH - Ratio of LF to HF

% frenquency ranges for powerband calculation
tot_range = [0.00, 3.99];
VLF_range = [0.00, 0.04];
LF_range = [0.05, 0.15];
HF_range = [0.16, 0.4];

tot_band_power = zeros(1,length(classes));
VLF_band_power = zeros(1,length(classes));
LF_band_power = zeros(1,length(classes));
HF_band_power = zeros(1,length(classes));
ratioLH = zeros(1,length(classes));
pVLF = zeros(1,length(classes));
pLF = zeros(1,length(classes));
pHF = zeros(1,length(classes));
tags = {'tot_bandpower' 'VLF_bandpower' 'LF_bandpower' 'HF_bandpower' 'ratioLH'};

for i = 1:5
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

data_out.frequency = table(tot_band_power', VLF_band_power', LF_band_power', ...
    HF_band_power', ratioLH', classes);

% Non-linear parameters calculated here are:
% -AppEn = approximate entropy
% -SaEn = sample entropy

dim = 2; % embedded dimension
tau = 0.1; %delay time lag
AppEn = zeros(1,length(classes));
SaEn = zeros(1,length(classes));
Alpha1 = zeros(1,length(classes));
Alpha2 = zeros(1,length(classes));
tags = {'App_En' 'Sa_En' 'alpha_1' 'alpha_2'};

for i = 1:5
    toll = .2 * NN_std(1,i); %tollerance
    AppEn(1,i) = ApEn( dim, toll, HRV_resample{i}); % Approximate entropy
    SaEn(1,i) = SampEn( dim, toll, HRV_resample{i}); %Sample Entropy
    %    [D2, Cm, epsilon] = corrdim(HRV_resample{i},dim,tau,epsilon,sloperange);
    %       m - embedding dimension
    %       tau - delay time lag
    [Alpha1(1,i), Alpha2(1,i)]=DFA_main(HRV_resample{i}); % Detrended fluctuation analysis
end

data_out.nonlinear = table (AppEn', SaEn', Alpha1', Alpha2', classes);
end