function [ results ] = process_timefreq_eeg( params, data, results )
% PROCESS_TIMEFREQ_EEG Performs time-frequency processing on EEG data for a
%   specific subject
%
% data:     Output from the preprocessing_eeg_3.m script
% results:  Struct to contain the results; if passed as an argument,
%           results will be appended to the struct under the field
%           "eeg.timefreq"
%

if nargin < 3
   results = []; 
end

subject = data.subject;




end

