%% Pre-process EEG data step 3 - Artifact rejection
%
% This script implements a processing pipeline for EEG data. It identifies
% ocular and zvalue artifacts, removes these, and interpolates.
%
% It requires step 2 to have been completed ("preprocessing_2_eeg.m"), as
% well as the eye preprocessing and processing pipelines 
% ("preprocessing_eye.m" and "processing_eye.m").
%
% Note: you need to have loaded a variable named "params", with all the
% required fields, for this pipeline to run properly. Use the
% "default_params_X.m" scripts to set your parameters, then run
% this script.
%
% Steps:
% 1. Load results of step 2 (ICA component removal)
% 2. Synchronise time series
% 3. Identify saccade intervals as artifacts
% 4. Identify z-score outliers as artifacts
% 5. Remove artifacts and save result
%

