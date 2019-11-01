function [ T, hdr ] = import_log( csv_file, fmt )

fid = fopen(csv_file);

hdr = fgetl(fid);
hdr = strrep(hdr,'"','');
hdr = strsplit(strtrim(hdr),',');
    
% User defined format 
T = textscan(fid, fmt, 'delimiter', ',');

fclose(fid);

end

