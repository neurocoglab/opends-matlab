function [ T ] = import_log_sim( csv_file, fmt )

fid = fopen(csv_file);

hdr = fgetl(fid);
hdr = strrep(hdr,'"','');
hdr = strrep(hdr,':','_');
hdr = strsplit(strtrim(hdr),',');

% User defined format 
T = textscan(fid, fmt, 'delimiter', ',');
C = [];

for i = 1 : length(T)
   E = T{i};
   if ~iscell(E)
       E = num2cell(E);
   end
   C = [C,E];
end

T = cell2table(C, 'VariableNames', hdr);

fclose(fid);

end

