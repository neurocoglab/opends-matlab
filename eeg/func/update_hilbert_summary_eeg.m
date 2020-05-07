function [ summary ] = update_hilbert_summary_eeg( results, summary )
%   UPDATE_HILBERT_SUMMARY_EEG Updates the summary variable for Hilbert
%       transforms; i.e., the summary over all subjects, which will be used to
%       perform statistical analysis

if isempty(summary)
    summary.subjects = {};
    summary.bands = results.eeg.hilbert.bands;
    summary.Fs = results.eeg.hilbert.Fs;
    summary.channels = {};
    summary.envelopes = {};
    summary.time = {}; 
    
end

summary.subjects = [summary.subjects {results.subject}];
summary.channels = [summary.channels {results.eeg.hilbert.channels}];
summary.envelopes = [summary.envelopes {results.eeg.hilbert.envelopes}];
summary.time = [summary.time {results.eeg.hilbert.time}];

end

