%% Main program, with classification and analysis.
%-----------------------Loading ECG Data----------------------------------%
% Version: Mihai, Vegard, Romek
% INCLUDE DESCRIPTION BEFORE HANDING IN
%
%
%
%-------------------------------------------------------------------------%

%prepare memory
close all; clear all; clc;

%% Load Data - ECG and time segments 

data_load;
%%
for subjectID = 1:2
   ecg=ecg_data(subjectID,~isnan(ecg_data(subjectID,:)));
   out(subjectID) = hrv_analysis (ecg, fs, startPoint(subjectID,:), classes);     
end

    




