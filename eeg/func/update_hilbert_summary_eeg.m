function [ summary ] = update_hilbert_summary_eeg( params, results, summary )
%   UPDATE_HILBERT_SUMMARY_EEG Updates the summary variable for Hilbert
%       transforms; i.e., the summary over all subjects, which will be used to
%       perform statistical analysis

if isempty(summary)
    summary.subjects = {};
    summary.bands = results.eeg.hilbert.bands;
    summary.Fs = results.eeg.hilbert.Fs;
    summary.N_cycles = [];
    summary.epochs.index.baseline = {};
    summary.epochs.index.overtake = {};
    summary.epochs.stats.baseline.mean = {};
    summary.epochs.stats.baseline.cycles.mean = {};
    summary.epochs.stats.overtake.mean = {};
    summary.epochs.stats.overtake.cycles.mean = {};
    summary.epochs.stats.cycles.mean = {};
    if params.sim.epochs.difficulty.apply
        summary.epochs.stats.overtake.difficulty.labels = params.sim.epochs.difficulty.labels;
        summary.epochs.stats.overtake.difficulty.mean = zeros(0,2);
        summary.epochs.overtake.index.difficulty = {};
    end
    if params.sim.epochs.outcomes.apply
        summary.epochs.stats.overtake.outcomes.labels = params.sim.epochs.outcomes.labels;
        summary.epochs.stats.overtake.outcomes.mean = zeros(0,2);
        summary.epochs.overtake.index.outcomes = {};
    end
    summary.channels = results.eeg.hilbert.channels;
end

summary.subjects = [summary.subjects {results.subject}];

summary.N_cycles = [summary.N_cycles;results.eeg.hilbert.N_cycles];

summary.epochs.index.baseline = [summary.epochs.index.baseline {results.eeg.hilbert.epochs.baselines}];
summary.epochs.index.overtake = [summary.epochs.index.overtake {results.eeg.hilbert.epochs.overtakes}];
summary.epochs.stats.baseline.mean = [summary.epochs.stats.baseline.mean {results.eeg.hilbert.epochs.stats.baseline.mean}];
summary.epochs.stats.overtake.mean = [summary.epochs.stats.overtake.mean {results.eeg.hilbert.epochs.stats.overtake.mean}];
summary.epochs.stats.cycles.mean = [summary.epochs.stats.cycles.mean {results.eeg.hilbert.epochs.stats.cycles.mean}];
summary.epochs.stats.baseline.cycles.mean = [summary.epochs.stats.baseline.cycles.mean ...
                                                {results.eeg.hilbert.epochs.stats.baseline.cycles.mean}];
summary.epochs.stats.overtake.cycles.mean = [summary.epochs.stats.overtake.cycles.mean ...
                                                {results.eeg.hilbert.epochs.stats.overtake.cycles.mean}];

if params.sim.epochs.difficulty.apply
    summary.epochs.overtake.index.difficulty = results.eeg.hilbert.epochs.overtake.difficulty;
    summary.epochs.stats.overtake.difficulty.mean = [summary.epochs.stats.overtake.difficulty.mean {results.eeg.hilbert.epochs.stats.overtake.difficulty.mean}];
end
if params.sim.epochs.outcomes.apply
    summary.epochs.overtake.index.outcome = results.eeg.hilbert.epochs.overtake.outcome;
    summary.epochs.stats.overtake.outcomes.mean = [summary.epochs.stats.overtake.outcomes.mean {results.eeg.hilbert.epochs.stats.overtake.outcomes.mean}];
end

end

