function [ ok ] = convert_log_sim( params, subject )
%CONVERT_LOG_SIMLOG Converts an xml-format simulation log to csv files

ok = true;

subj_dir = sprintf('%s/%s/%s', params.io.input_dir, params.sim.sub_dir, subject);
log_file = sprintf('%s/%s%s.log', subj_dir, params.sim.convert.prefix, subject);
exec = fullfile(params.general.matlab_dir, 'bin', params.sim.convert.exec);
outdir = sprintf('%s/%s/%s', params.io.output_dir, subject, params.sim.sub_dir);

if ~exist(outdir, 'dir')
   mkdir(outdir); 
end

a = strfind(exec,'/');
a = a(end);
cdir = exec(1:a);
cfile = exec(a+1:end);
if ~ispc, cfile = ['./' cfile]; end

cmd = sprintf('cd "%s"; %s "%s" assets/opends/%s.xml "%s" -addserial', cdir, ...
                                            cfile, ...
                                            log_file, ...
                                            params.sim.convert.filter, ...
                                            outdir);

[status,result] = system(cmd);
if status ~= 0 || contains(result,'Exception')
    warning('Error converting simlog for %s: %s', subject, result);
    ok = false;
else
%     fid = fopen(flag_file,'w');
%     fclose(fid);
%     fprintf('Done.\n');
end    


end

