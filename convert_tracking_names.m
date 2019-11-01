track_dir = '/Volumes/AndrewElements/data/driving/tracking';

files = dir(track_dir);
files = files([files.isdir]);
files = files(3:end);

for i = 1 : length(files)
    subject = files(i).name;
    zipfile = sprintf('%s/%s/%s-tracking.zip', track_dir, subject, subject);
    if exist(zipfile, 'file')
       newfile = sprintf('%s/%s/tracking-%s.zip', track_dir, subject, subject);
       movefile(zipfile, newfile)
       fprintf('Changed: %s\n', zipfile)
    end
end