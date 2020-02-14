function [ ok ] = convert_eyelog_smi( subject, params )
%convert_eyelog_smi Converts the ASCII version of an SMI eye tracker log 
%                   (or multiple logs) to CSV files
%

ok = true;

subj_dir = sprintf('%s/%s/%s', params.io.input_dir, params.eye.sub_dir, subject);
output_dir = sprintf('%s/%s/%s', params.io.output_dir, subject, params.eye.sub_dir);
temp_dir = sprintf('%s/tmp', output_dir);

if ~exist(output_dir, 'dir')
   mkdir(output_dir); 
end

if ~exist(temp_dir, 'dir')
   mkdir(temp_dir); 
end
                               
log_file = sprintf('%s/%s_samples.txt', subj_dir, subject);

log_files = [];
parts = [];

if ~exist(log_file, 'file')
    [a,fn,ext] = fileparts(log_file);
    log_file = sprintf('%s/%s_part1%s', a, fn, ext);
    k = 2;
    while exist(log_file, 'file')
        log_files = [log_files {log_file}];
        log_file = sprintf('%s/%s_part%d%s', a, fn, k, ext);
        k = k + 1;
    end
else
    log_files = {log_file};
end

if isempty(log_files)

  log_file = sprintf('%s/%s_samples.txt', temp_dir, ...
                                          subject);

  zip_file = sprintf('%s/%s%s.zip', subj_dir, ...
                                    params.eye.convert.prefix, ...
                                    subject);

  try
      if exist(zip_file, 'file')
         % Extract samples log from zip file
         unzip(zip_file, temp_dir);

         % Is is in the current directory?
         lfile = dir(sprintf('%s/*Samples.txt', temp_dir));

         if ~isempty(lfile)
             lfile = lfile(1);
             movefile(sprintf('%s/%s', temp_dir, lfile.name), ...
                      sprintf('%s/tmp', temp_dir));
             delete(sprintf('%s/*.txt', temp_dir));
             movefile(sprintf('%s/tmp', temp_dir), ...
                      log_file);
             log_files = {log_file};
         else
             % Otherwise, find subdirectory containing subject
             % ID (if this doesn't exist we have a problem)
             subdirs = dir(temp_dir);
             subdirs = subdirs(3:end);
             subdirs = subdirs([subdirs.isdir]);
             subdir = [];
             for s = 1 : length(subdirs)
                if strfind(subdirs(s).name, subject)
                   subdir = subdirs(s).name;
                end
             end
             if ~isempty(subdir)
                 lfiles = dir(sprintf('%s/%s/*Samples.txt', temp_dir, subdir));

                 for f = 1 : length(lfiles)
                     lfile = lfiles(f);
                     fname = lfile.name;
                     idx = strfind(fname,'part');
                     if idx > 0
                         k = fname(idx+4);
                         lfile2 = sprintf('%s/%s/%s', temp_dir, subdir, fname);
                         out_file = log_file;
                         if length(lfiles) > 1
                             [a,fn,ext] = fileparts(out_file);
                             out_file = sprintf('%s/%s_part%s%s', a, fn, k, ext);
                         end
                         movefile(lfile2, out_file);
                         log_files = [log_files {out_file}];
                         parts = [parts {k}];
%                                      fprintf('Created %s\n', out_file);
                     else
                        warning('No part in sample file name: %s', fname)
                     end
                 end
                 rmdir(sprintf('%s/%s/%s', temp_dir, subject),'s');
             end
         end
      end
  catch err
      % This will get handled by the next statement 
      a=0;
  end
end

if isempty(log_files)
   fprintf('\nLog or zip file not found for %s; skipping subject.\n', subject);
   ok = false;
else
   prefix = params.eye.convert.prefix;
   for j = 1 : length(log_files)
       log_file = log_files{j};
       if length(log_files) > 1
           idx = strfind(log_file,'part');
           if idx > 0
               k = log_file(idx+4);
           end
           prefix = sprintf('%spart%s', params.eye.convert.prefix, k);
       end
       
       exec = fullfile(params.general.matlab_dir, 'bin', params.eye.convert.exec_smi);

       a = strfind(exec,'/');
       a = a(end);
       cdir = exec(1:a);
       cfile = exec(a+1:end);
       if ~ispc, cfile = ['./' cfile]; end
       cmd = sprintf('cd "%s"; pwd; %s "%s" -o "%s" -columns "%s" -prefix %s -clobber', ...
                        cdir, ...
                        cfile, ...
                        log_file, ...
                        output_dir, ...
                        params.eye.convert.columns, ...
                        prefix);
%                     fprintf('%s\n',cmd);
        [status,result] = system(cmd);
        if status ~= 0
            fprintf('\nError converting log for %s.\n%s', subject, result);
            ok = false;
        end
   end
   
   if ok
      % Clean up
      rmdir(temp_dir, 's');
   else
       
%       fprintf('Errors encountered.\n');
   end
end
           

end

