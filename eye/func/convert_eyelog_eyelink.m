function [ ok ] = convert_eyelog_eyelink( subject, params )
%convert_eyelog_eyelink Converts the ASCII version of an Eyelink eye tracker log 
%                       (or multiple logs) to CSV files
%



rdir = params.root_dir;
if ~isempty(params.data_dir)
  rdir = sprintf('%s/%s', rdir, params.data_dir);
end

input_dir = sprintf('%s/%s/%s', rdir, ...
                                   params.convert.input_dir, ...
                                   subject);

output_dir = sprintf('%s/%s/%s', rdir, ...
                                   params.output_dir, ...
                                   subject);
 
if ~exist(output_dir, 'dir')
   mkdir(output_dir); 
end

if ~isempty(params.convert.edf2asc)
   fprintf(' Converting from EDF to ASC...\n');
   edf_files = dir(sprintf('%s/*.edf', input_dir));
   
   for i = 1 : length(edf_files)
      filei = sprintf('%s/%s', edf_files(i).folder, edf_files(i).name);
      fileasc = [filei(1:end-4) '.asc'];
      if exist(fileasc, 'file')
          fprintf(' already converted %s\n', edf_files(i).name);
      else
          cmd = sprintf('"%s" -p "%s" "%s"', params.convert.edf2asc, edf_files(i).folder, filei);
          [status, message] = system(cmd);
          if status ~= 255
              fprintf(' Problem converting %s: %s\n', filei, message);
              ok = false;
              return;
          end
      end
   end
    
end

try
                               
    log_file = sprintf('%s/%s%s.asc', input_dir, params.convert.prefix, subject);
    if ~exist(log_file, 'file')
        zip_file = sprintf('%s/ods%s.zip', input_dir, subject);

        if ~exist(zip_file, 'file')
            warning('No data found for %s! Skipping.', subject);
            ok = false;
            return;
        end

        unzip(zip_file, input_dir);
        log_file = sprintf('%s/ods%s.asc', input_dir, subject);

        if ~exist(log_file, 'file')
            warning('No data found in zip file for %s! Skipping.', subject);
            ok = false;
            return;
        end

    end

catch err
    warning('Exception: %s', err);
    ok = false;
    return;
end

% Read Eyelink log and convert to CSV
cfg.dataset = log_file;
[~,data_eye] = evalc('ft_preprocessing(cfg);');

samples_out = sprintf('%s/%s_samples.csv', output_dir, subject);

[fid_out, message] = fopen(samples_out, 'w+');
if fid_out < 0
    warning('Could not open CSV output file %s, with error: %s', samples_out, message);
    fclose all;
    ok = false;
    return;
end

% Write samples
fprintf(fid_out, 'Time,Type,GazeX,GazeY,PupilDiam\n');
T = data_eye.trial{1};
for i = 1 : size(data_eye.trial{1},2)
    fprintf(fid_out, '%d,%s,%1.4f,%1.4f,%1.4f\n', T(1,i), 'SMP', T(2,i), T(3,i), T(4,i)); 
end

clear T data_eye;
fclose(fid_out);

[~,event_eye] = evalc('ft_read_header(cfg.dataset);');

messages_out = sprintf('%s/%s_messages.csv', output_dir, subject);
[fid_out, message] = fopen(messages_out, 'w+');
if fid_out < 0
    warning('Could not open CSV output file %s, with error: %s', messages_out, message);
    fclose all;
    ok = false;
    return;
end

% Write messages
fprintf(fid_out, 'Time,Trial,LogId,EventType\n');
for i = 1 : length(event_eye.orig.msg)
    % MSG	307833 EVENT[id=2;type=SimulatorEnded]
   msg = event_eye.orig.msg{i};
   if contains(msg, 'EVENT[')
       parts = strsplit(msg,[{'\t'},{' '}]);
       t = parts{2};
       parts = strsplit(parts{3}(7:end-1),[{'='},{';'}]);
       fprintf(fid_out, '%s,1,%s,%s\n', t, parts{2}, parts{4});
   end
    
end

fclose(fid_out);

ok = true;

end

