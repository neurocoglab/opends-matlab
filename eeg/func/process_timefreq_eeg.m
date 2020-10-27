function [ results ] = process_timefreq_eeg( params, data, results )
% PROCESS_TIMEFREQ_EEG Performs time-frequency processing on EEG data for a
%   specific subject
%
% Arguments:
% data:     Output from the preprocessing_eeg_3.m script
% results:  Output of "processing_eye.m"
%
% Returns:
% results:  a struct containing time
%

if nargin < 3
   results = []; 
end

subject = data.subject;
outdir = sprintf( '%s/%s', params.io.output_dir, subject );
timefreq_dir = sprintf( '%s/timefreq', outdir );
if ~exist( timefreq_dir, 'dir' )
   mkdir( timefreq_dir );
end

results = get_timefreq_trials_eeg( params, data, results );

% outdir = sprintf( '%s/%s', params.io.output_dir, subject );
% fft_file = sprintf('%s/results_fft_eeg.m', outdir);
% 
% if exist(fft_file, 'file') && ~params.eeg.timefreq.regen_fft
%     Load existing FFT results
%     load( fft_file );
%     
% else
%     
%     Compute time/frequency for each round
%     cfg = [];
%     cfg.toi = [];
%     cfg.Zscore = false;
%     cfg.dt = params.eeg.timefreq.dt;
%     cfg.foi = params.eeg.timefreq.foi;
%     cfg.wavelet = params.eeg.timefreq.wavelet;
% 
%     results.eeg.timefreq = process_wavelet_eeg(cfg, results.eeg.timefreq, [], params.general.debug);
% 
%     Save FFT results
%     if params.general.debug
%        fprintf(' * Saving FFT to %s...', fft_file);
%     end
%     save( fft_file, 'results' );
%     if params.general.debug
%        fprintf('done.\n');
%     end
% end
% 
% Perform statistical analyses
% results = analyze_timefreq_eeg( params, data, results );
% 
% 

end

