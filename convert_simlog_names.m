simlog_dir = '/Volumes/AndrewElements/data/driving/driving';

files = dir(simlog_dir);
files = files([files.isdir]);
files = files(3:end);

for i = 1 : length(files)
    subject = files(i).name;
    logfile = sprintf('%s/%s/%s-simlog.log', simlog_dir, subject, subject);
    if exist(logfile, 'file')
       newfile = sprintf('%s/%s/simlog-%s.log', simlog_dir, subject, subject);
       movefile(logfile, newfile)
       fprintf('Changed: %s\n', logfile)
    end
end