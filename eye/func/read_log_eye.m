function [ data ] = read_log_eye( params, subject )
%read_log_eye Reads converted eye data and saves in Matlab format

subj_dir = sprintf('%s/%s/%s', params.io.output_dir, subject, params.eye.sub_dir);
prefix = [params.eye.convert.prefix subject];
input_files = [];
input_file = sprintf('%s/%s_samples.csv', subj_dir, prefix);

if exist(input_file, 'file')
   input_files = {input_file};
else
   input_file = sprintf('%s/%s-part1_samples.csv', subj_dir, prefix);
   k = 2;
   while exist(input_file, 'file')
       input_files = [input_files {input_file}];
       input_file = sprintf('%s/%s-part%d_samples.csv', subj_dir, prefix, k);
       k = k + 1;
   end
end

Fs = NaN;
data = []; data_i = [];
data.subject = subject;
data_i.subject = subject;

for j = 1 : length(input_files)

    input_file = input_files{j};
    T = import_samples_eye(input_file);

    data_i.eye.t = T{1}(2:end); % cell2mat(T{1}(2:end));
    data_i.eye.pos_x = T{3}(2:end); % cell2mat(T(2:end,3));
    data_i.eye.pos_y = T{4}(2:end); % cell2mat(T(2:end,4));
    data_i.eye.diam = T{5}(2:end); % cell2mat(T(2:end,5));

    clear T;

    % Convert double to single
    data_i.eye.pos_x = single(data_i.eye.pos_x);
    data_i.eye.pos_y = single(data_i.eye.pos_y);
    data_i.eye.diam = single(data_i.eye.diam);

    % Make data compatible
    if strcmp(params.eye.convert.format, 'smi')
        data_i.eye.t = data_i.eye.t / 1000;
    end
    if strcmp(params.eye.convert.format, 'eyelink')
        data_i.eye.diam = data_i.eye.diam / 10;
    end
    
    % Record start time
    if j == 1
        data.eye.t_start = data_i.eye.t(1);
        data_i.eye.t_start = data.eye.t_start;
    else
        data_i.eye.t_start = data.eye.t_start + data_i.eye.t(1);
    end

%     data_i.eye.t = data_i.eye.t - data.eye.t_start;
    data_i.eye.t = single(data_i.eye.t);
    data_i.eye.Fs = Fs;

    if j == 1
       data = data_i;
    else
       % Append to existing data
       data.eye.t = [data.eye.t;data_i.eye.t];
       data.eye.pos_x = [data.eye.pos_x;data_i.eye.pos_x];
       data.eye.pos_y = [data.eye.pos_y;data_i.eye.pos_y];
       data.eye.diam = [data.eye.diam;data_i.eye.diam];
    end

end

data.eye.Fs = params.eye.Fs;

% Find gaps in sampling
if params.eye.gaps.apply
    data = remove_gaps_eye( data, params );
end

% Load messages from eye tracker
input_file = sprintf('%s/%s_messages.csv', subj_dir, prefix);

if exist(input_file, 'file')
   [messages, hdr] = import_messages_eye(input_file);
else
   % There are probably multiple parts; need to glue these together
   input_file = sprintf('%s/%s-part1_messages.csv', subj_dir, prefix);
   if ~exist(input_file, 'file')
       error('No message log found for subject %s', subject);
   end
   
   hdr = [];
   messages = [];
   k = 1;
   while exist(input_file, 'file')
       [messages_i, hdr] = import_messages_eye(input_file);
       if k == 1
           messages = messages_i;
       else
           for j = 1 : 4
               messages(j) = {[messages{j};messages_i{j}]};
           end
       end
       k = k + 1;
       input_file = sprintf('%s/%spart%d_messages.csv', subj_dir, prefix, k);
       
   end
   
end

messages = [num2cell(messages{:,1}),num2cell(messages{:,2}), ...
               num2cell(messages{:,3}),num2cell(messages{:,4})];
messages = cell2table(messages, 'VariableNames', hdr);

if strcmp(params.eye.convert.format, 'smi')
    % Convert SMI (microseconds) to milliseconds
    messages.Time = messages.Time / 1000;
end

data.eye.log.messages = messages;


end

