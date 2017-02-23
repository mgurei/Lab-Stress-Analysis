%-----------------------Loading ECG Data----------------------------------%
% Version: George
%-- This file loads the ECG data and the starting points of every
% stage for each subject (see paper on Acute Stress by Tanev et all)

%% ************** Load Data
ecg_data=nan(10,max(max(startPoint))+fs*60*15);
for i=1:10
    data=my_edfRead([names{i} '.edf']);  
    ecg_data(i,1:length(data(1,:)))=data(1,:);
end
save('s_ecg_10_subjects','ecg_data');

%% Start sample point of the four stages (columns) for every subject (rows)
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
%% Example Detect QRS
load s_b_coeff.mat;
load s_ecg_10_subjects.mat;
ecg=ecg_data(1,~isnan(ecg_data(1,:)));
qrs=s_detect_qrs(ecg,b_low,b_high,b_deriv,b_avg,delay);
figure;
plot(1:length(qrs),5000*qrs,1:length(ecg),ecg);

%% Extracting HRV Parameters

%-- HRV Time Domain Measures
i=1;
meanNN(i)=mean(HRV); % Mean NN times time
sdNN(i)=std(HRV); % Standard deviation of NN times
% RMS of differences between adjacent NN intervals
RMSSD(i)=sqrt(mean(diff(HRV).^2));
SDSD(i)=std(diff(HRV)); % STD of differences between adjacent NN intervals
NN50(i)=sum(0.05<diff(HRV));
pNN50(i)=NN50(i)/length(HRV)*100;
%-- HRV Frequency Measures
% Total HRV Positive frequency range
Total=find(freq_axis>=0.00 & freq_axis<f_resample/2);
Total_sum(i)=sum(abs(HRV_psd(Total)));
% VLF - HRV from 0.00 to 0.04
VLF=find(freq_axis>=0.00 & freq_axis<0.05);
VLF_sum(i)=sum(abs(HRV_psd(VLF)));
% LF - HRV from 0.05 to 0.15 Hz
LF=find(freq_axis>=0.05 & freq_axis<0.15);
LF_sum(i)=sum(abs(HRV_psd(LF)));
LF_norm(i)=LF_sum(i)/(Total_sum(i)-VLF_sum(i));
% HF - HRV from 0.16 to 0.40 Hz
HF=find(freq_axis>=0.15 & freq_axis<=0.4);
HF_sum(i)=sum(abs(HRV_psd(HF)));
HF_norm(i)=HF_sum(i)/(Total_sum(i)-VLF_sum(i));
% Ratio of LF to HF
LF_HF(i)=LF_norm(i)/HF_norm(i);

%-- Non-Linear Measures
% r value necessary for entropy calculations (recommended by Chon et
% al.)
rch=(-0.036+0.26*sqrt(SDSD(i)/sdNN(i)))/nthroot(length(HRV)/1000,4);
%rch_re=(-0.036+0.26*sqrt(std(diff(HRV_resample))/std(HRV_resample)))/nthroot(length(HRV_resample),4);
% Approximate Entropy, r=0.2sdNN
AppEn02(i)=ApEn(2,0.2*sdNN(i),HRV_resample,1); % sdNN is of HRV - NOT HRV resample
% Approximate Entropy, r=rch
AppEnCho(i)=ApEn(2,rch,HRV_resample,1); % sdNN is of HRV - NOT HRV resample
SampleEn02(i)=SampEn(2,0.2*sdNN(i),HRV,1); % Sample Entropy
SampleEnCho(i)=SampEn(2,rch,HRV,1);
[Alph1(i) Alph2(i)]=DFA_main(HRV); % Dentreded Fluctuation Analysis
SD1(i)=sqrt(SDSD(i)^2/2); % Poincarre Plot SD1
SD2(i)=sqrt(2*sdNN(i)^2-SDSD(i)^2/2); %  Poincarre Plot SD2
% Correlation Dimension
m=10; tau=1; % Embedding Dimension and Delay respectively
eps=0.001:0.001:2; % Radius threshold vector differences in embedded dimension 
sloperange=[-8 -5]; % Range of log2(Cm) for slope calculation
[D2(i), Cm, epsilon]=corrdim(HRV,m(end),tau,eps,sloperange);
