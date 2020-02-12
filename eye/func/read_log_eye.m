function [ data ] = read_eye_data( subject, params )
%read_eye_data Reads convert eye data and saves in Matlab format


input_files = [];
input_file = sprintf('%s/%s/%s/%s/%s_samples.csv', params.root_dir, params.data_dir, params.output_dir, subject, subject);

if exist(input_file, 'file')
   input_files = {input_file};
else
   input_file = sprintf('%s/%s/%s/%s/%s_part1_samples.csv', params.root_dir, params.data_dir, params.output_dir, subject, subject);
   k = 2;
   while exist(input_file, 'file')
       input_files = [input_files {input_file}];
       input_file = sprintf('%s/%s/%s/%s_part%d_samples.csv', params.root_dir, params.data_dir, params.output_dir, subject, subject, k);
       k = k + 1;
   end
end

Fs = NaN;
data = []; data_i = [];
data.subject = subject;
data_i.subject = subject;

for j = 1 : length(input_files)

    input_file = input_files{j};
    T = import_samples(input_file);

    data_i.eye.t = T{1}(2:end); % cell2mat(T{1}(2:end));
    data_i.eye.pos_left_x = T{3}(2:end); % cell2mat(T(2:end,3));
    data_i.eye.pos_left_y = T{4}(2:end); % cell2mat(T(2:end,4));
    data_i.eye.diam_left = T{5}(2:end); % cell2mat(T(2:end,5));

    clear T;

    % Convert double to single
    data_i.eye.pos_left_x = single(data_i.eye.pos_left_x);
    data_i.eye.pos_left_y = single(data_i.eye.pos_left_y);
    data_i.eye.diam_left = single(data_i.eye.diam_left);

    % Make data compatible
    if strcmp(params.tracker_type, 'smi')
        data_i.eye.t = data_i.eye.t / 1000;
    end
    if strcmp(params.tracker_type, 'eyelink')
        data_i.eye.diam_left = data_i.eye.diam_left / 10;
    end
    
    % Convert time to ms relative to start
    if j == 1
        data.eye.t_start = data_i.eye.t(1);
        data_i.eye.t_start = data.eye.t_start;
    else
        data_i.eye.t_start = data.eye.t_start + data_i.eye.t(1);
    end

    data_i.eye.t = data_i.eye.t - data.eye.t_start;
    data_i.eye.t = single(data_i.eye.t);
    data_i.eye.Fs = Fs;

    if j == 1
       data = data_i;
    else
       % Append to existing data
       data.eye.t = [data.eye.t;data_i.eye.t];
       data.eye.pos_left_x = [data.eye.pos_left_x;data_i.eye.pos_left_x];
       data.eye.pos_left_y = [data.eye.pos_left_y;data_i.eye.pos_left_y];
       data.eye.diam_left = [data.eye.diam_left;data_i.eye.diam_left];
    end

end

data.eye.Fs = params.Fs;

% Find gaps in sampling
if params.remove_gaps
    data.eye = remove_gaps( data.eye, params.blink );
end

end

