function [ score ] = get_subject_score( params, subject )

% Open subject simlog
subj_dir = sprintf('%s/%s/%s', params.io.input_dir, params.sim.sub_dir, subject);
log_file = sprintf('%s/%s%s.log', subj_dir, params.sim.convert.prefix, subject);
if ~exist(log_file, 'file')
   score = nan;
   return;
end

% Get final score
fid = fopen( log_file );

score = nan;
tline = fgetl(fid);

while ischar(tline)
    if contains(tline, '<Game:AccumulatedScore>')
       i0 = strfind(tline, '>');
       i1 = strfind(tline, '<');
       if ~(isempty(i0) || length(i1) < 2)
           score = str2double(tline(i0(1)+1:i1(2)-1));
       end
    end
    
    tline = fgetl(fid);
end

end

