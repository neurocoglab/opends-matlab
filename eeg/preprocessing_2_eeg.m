%% Pre-process EEG data step 2 - Interactive ICA component rejection
%
% This script allows the user to interactively select and reject individual
% ICA components. It requires step 1 to have been completed
% ("preprocessing_1_eeg.m").
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
% 4. Remove components and save result as "eeg-{uid}_preproc2.mat"
%
%