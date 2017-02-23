%% Main program, with classification and analysis.
%-----------------------Loading ECG Data----------------------------------%
% Version: Mihai, Vegard, Romek
% INCLUDE DESCRIPTION BEFORE HANDING IN
%-------------------------------------------------------------------------%

%tidying up
close all; clear all; clc;

%% Load Data - ECG and time segments

data_load;
%%

%Initialization of feature vectors.
feature_names = {'NN_mean'; 'NN_std'; 'NN_diff_rms'; 'NN_diff_std'; 'NN50'; 'pNN50'; ...
    'tot_band_power'; 'VLF_band_power'; 'LF_band_power'; 'HF_band_power'; 'ratioLH'; ...
    'AppEn'; 'SaEn'; 'Alpha1'; 'Alpha2'};

NN_mean = 0;
NN_std = 0;
NN_diff_rms = 0;
NN_diff_std = 0;
NN50 = 0;
pNN50 = 0;
tot_band_power = 0;
VLF_band_power = 0;
LF_band_power = 0;
HF_band_power = 0;
ratioLH = 0;
AppEn = 0;
SaEn = 0;
Alpha1 = 0;
Alpha2 = 0;

class_vector = [];
out = struct;

for subjectID = 1:10
    ecg=ecg_data(subjectID,~isnan(ecg_data(subjectID,:)));
    out = hrv_analysis (ecg, fs, startPoint(subjectID,:), classes);
    
    % Reorganizing temporal features
    NN_mean = [NN_mean; out.time.Var1];
    NN_std = [NN_std; out.time.Var2];
    NN_diff_rms = [NN_diff_rms; out.time.Var3];
    NN_diff_std = [NN_diff_std; out.time.Var4];
    NN50 = [NN50; out.time.Var5];
    pNN50 = [pNN50; out.time.Var6];
    
    % Reorganizing frequency features
    tot_band_power = [tot_band_power; out.frequency.Var1];
    VLF_band_power = [VLF_band_power; out.frequency.Var2];
    LF_band_power = [LF_band_power; out.frequency.Var3];
    HF_band_power = [HF_band_power; out.frequency.Var4];
    ratioLH = [ratioLH; out.frequency.Var5];
    
    % Reorganizing non-linear features
    AppEn = [AppEn; out.nonlinear.Var1];
    SaEn = [SaEn; out.nonlinear.Var2];
    Alpha1 = [Alpha1; out.nonlinear.Var3];
    Alpha2 = [Alpha2; out.nonlinear.Var4];
    
    class_vector = [class_vector; classes];
end

NN_mean(1) = [];
NN_std(1) = [];
NN_diff_rms(1) = [];
NN_diff_std(1) = [];
NN50(1) = [];
pNN50(1) = [];
tot_band_power(1) = [];
VLF_band_power(1) = [];
LF_band_power(1) = [];
HF_band_power(1) = [];
ratioLH(1) = [];
AppEn(1) = [];
SaEn(1) = [];
Alpha1(1) = [];
Alpha2(1) = [];
%%
feature_matrix = [NN_mean, NN_std, NN_diff_rms, NN_diff_std, NN50, pNN50, ...
    tot_band_power, VLF_band_power, LF_band_power, HF_band_power, ratioLH, ...
    AppEn, SaEn, Alpha1, Alpha2];

save ('feature_matrix.mat', 'feature_matrix')
save ('class_vector.mat', 'class_vector')


%%
load('feature_matrix.mat')
load('class_vector.mat')

cp = cvpartition(class_vector, 'KFold')

%
bayes_fit = fitcnb(feature_matrix(:,1:15), class_vector, 'DistributionNames',...
    'kernel', 'Kernel','normal');
fit_Class = resubPredict(bayes_fit);
[ldaResubCM,grpOrder] = confusionmat(class_vector, fit_Class)
fit_ResubErr = resubLoss(bayes_fit)
fit_CV = crossval(bayes_fit, 'CVPartition',cp);
fit_CVErr = kfoldLoss(fit_CV)






