function [ ok ] = convert_eyelog_tobii( subject, params )
%convert_eyelog_tobii Converts a CSV-format log generated from a Tobii eye tracker
%                       (using a custom Python server script)
%

subj_dir = sprintf('%s/%s/%s', params.io.original_dir, params.eye.sub_dir, subject);
output_dir = sprintf('%s/%s/%s', params.io.output_dir, subject, params.eye.sub_dir);

if exist(output_dir, 'dir')
   rmdir( output_dir, 's' );
end
mkdir(output_dir);

% Read Tobii log
log_file = sprintf('%s/%s.csv', subj_dir, subject);
opts = detectImportOptions(log_file);
opts.VariableTypes(find(strcmp(opts.VariableNames, 'Message'))) = {'char'};
opts.VariableNamingRule = 'preserve';

log_table = readtable(log_file, opts);

% Remove invalid rows
[log_table, rem_count] = clean_up_data( log_table );

% Sort log entries by simulation time
log_table = sortrows(log_table, 'SimulatorTime');

% Filter for samples
prefix = [params.eye.convert.prefix subject];
samples_out = sprintf('%s/%s_samples.csv', output_dir, prefix);

[fid_smp, message] = fopen(samples_out, 'w+');
if fid_smp < 0
    fprintf('Could not open CSV output file %s, with error: %s', samples_out, message);
    fclose all;
    ok = false;
    return;
end

messages_out = sprintf('%s/%s_messages.csv', output_dir, prefix);
[fid_msg, message] = fopen(messages_out, 'w+');
if fid_msg < 0
    fprintf('Could not open CSV output file %s, with error: %s', messages_out, message);
    fclose all;
    ok = false;
    return;
end

if rem_count > 0
    warning('Removed %d invalid/corrupt rows from Tobii log; see log file.\n', rem_count);
end

% Write samples and events
fprintf(fid_smp, 'Time,Type,GazeX,GazeY,PupilDiam\n');
fprintf(fid_msg, 'Time,Trial,LogId,EventType\n');
for i = 1 : height(log_table)
    if strcmp(log_table.SampleType{i}, 'Gaze')
        fprintf(fid_smp, '%d,%s,%1.4f,%1.4f,%1.4f\n', log_table.SimulatorTime(i), 'SMP', ...
                                                      mean([log_table.GazeLeftX(i) log_table.GazeRightX(i)]), ...
                                                      mean([log_table.GazeLeftY(i) log_table.GazeRightY(i)]), ...
                                                      mean([log_table.PupilLeft(i) log_table.PupilRight(i)]));
    elseif strcmp(log_table.SampleType{i}, 'Message')
        msg = log_table.Message{i};
        if startsWith(msg, 'EVENT ')
            % "EVENT id=6;type=SimulatorStarted;millis=1711535298990"
           parts = strsplit(msg,[{'\t'},{' '}]);
           parts = strsplit(parts{2},[{'='},{';'}]);
           fprintf(fid_msg, '%d,1,%s,%s\n', log_table.SimulatorTime(i), parts{2}, parts{4});
        end
    end
end

fclose(fid_smp);
fclose(fid_msg);

ok = true;

    

end

function [log_table, rem_count] = clean_up_data( log_table )
        % Clean up bad lines (likely due to synchronization issues)

        rem_count = 0;

        % 1. Extra variables?
        idx_rem = [];
        for k = 1 : 20
            varname =  sprintf('ExtraVar%d',k);
            if any(strcmp(log_table.Properties.VariableNames,varname))
                idx_rem = [idx_rem;find(~strcmp(log_table.(varname),''))];
            end
        end
        
        if ~isempty(idx_rem)
            idx_rem = unique(idx_rem);
            idx_keep = true(height(log_table),1);
            idx_keep(idx_rem) = false;
            log_table = log_table(idx_keep,:);
            rem_count = length(idx_rem);
        end

        % 2. Time is extreme outlier
        time_vars = {'SystemTime','SimulatorTime'};
        for k = 1 : length(time_vars)
            tvz = log_table.(time_vars{k});
            idx_rem = isnan(tvz);
            if sum(idx_rem) > 0
                log_table = log_table(~idx_rem,:);
                rem_count = rem_count + sum(idx_rem);
                tvz = log_table.(time_vars{k});
            end
            
            [~,idx_sort] = sort(tvz);
            D = diff(idx_sort);
            idx_rem = false(height(log_table),1);
            idx_rem(idx_sort(abs(D) > 10E3)) = true;
            if sum(idx_rem) > 0
                log_table = log_table(~idx_rem,:);
                rem_count = rem_count + sum(idx_rem);
            end
        end
        
        % 3. Not Gaze or Message
        idx_rem = ~(strcmp(log_table.SampleType,'Gaze') | strcmp(log_table.SampleType,'Message'));
        if sum(idx_rem) > 0
            log_table = log_table(~idx_rem,:);
            rem_count = rem_count + sum(idx_rem);
        end


    end
