function [ summary ] = analyze_hilbert_eeg( params, summary )
% ANALYZE_HILBERT_EEG Performs statistical analyses on Hilbert envelope
%   output
%   
% Note: Hilbert processing must already have been performed
%

summary = [];

% 1. Epochs
%    Does Hilbert envelope magnitude differ between epochs? 

for i = 1 : length(summary.subjects)
    
    subject = summary.subjects{i};
    
    %    1.1. Passing vs. baseline
    envelopes = summary.envelopes{i};

    %    1.2. Good vs. poor outcome


    
    
    
end


% 2. Rounds
%    Does Hilbert envelope differ across rounds?
%    2.1. Overall

%    2.2. Within epochs (passing/baseline)



% 3. Correlation with pupil?
%    Does Hilbert envelope covary with pupil diameter and 1st temporal 
%     derivative?
%    3.1. Full time series


%    3.2. Moving window (dynamic)


%    3.3. Within epochs (passing/baseline)






end

