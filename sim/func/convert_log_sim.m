function [ ok, score ] = convert_log_sim( params, subject )
%CONVERT_LOG_SIMLOG Converts an xml-format simulation log to csv files

ok = true;

subj_dir = sprintf('%s/%s/%s', params.io.original_dir, params.sim.sub_dir, subject);
log_file = sprintf('%s/%s%s.log', subj_dir, params.sim.convert.prefix, subject);
exec = fullfile(params.general.matlab_dir, 'bin', params.sim.convert.exec);
outdir = sprintf('%s/%s/%s', params.io.output_dir, subject, params.sim.sub_dir);

if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

exec = strrep(exec,'\','/');
a = strfind(exec,'/');
a = a(end);
exec = strrep(exec,'/',filesep);
cdir = exec(1:a);
cfile = exec(a+1:end);

if ismac || isunix
    cfile = ['./' cfile];
    cmd = sprintf('cd "%s"; %s "%s" assets%sopends%s%s.xml "%s" -addserial', cdir, ...
                                                cfile, ...
                                                log_file, ...
                                                filesep, filesep, ...
                                                params.sim.convert.filter, ...
                                                outdir);
else
    cmd = sprintf('cd "%s" & %s "%s" assets%sopends%s%s.xml "%s" -addserial', cdir, ...
                                                cfile, ...
                                                log_file, ...
                                                filesep, filesep, ...
                                                params.sim.convert.filter, ...
                                                outdir);
end

if params.general.debug
   fprintf('\nDEBUG: Running: %s\n', cmd); 
end
[status,result] = system(cmd);
if status ~= 0 || contains(result,'Exception')
    warning('Error converting simlog for %s: %s', subject, result);
    ok = false;
    return;
end

% Compile final scores for all subjects
score = get_subject_score( params, subject, log_file );

% Write to file
fid = fopen(sprintf('%s/final_score.csv', outdir),'w');
fprintf(fid, '%d', score);
fclose(fid);

end

