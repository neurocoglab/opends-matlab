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
% "set_preprocessing_params_eeg.m" script to set your parameters, then run
% this script.
%
% Steps:
% 1. Load ICA result
% 2. Show ICA results
% 3. Prompt for components to reject
% 4. Remove components and save result as "eeg-{uid}_preproc3.mat"
%

