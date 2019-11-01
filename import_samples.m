function [ T, hdr ] = import_samples( csv_file )

% [~, result] = system( ['wc -l ', csv_file] );
% result=strsplit(strtrim(result),' ');
% N_lines = str2num(result{1});
% 
% fid = fopen(csv_file);
% 
% hdr = fgetl(fid);
% hdr = strsplit(strtrim(hdr),',');
% N_cols = length(hdr);
% 
% T = cell(N_lines, N_cols);
% 
% for i = 1 : N_lines
%     
%     %800181281,"SMP",519.0766,428.2381,21.528,21.528
%     T(i,:) = textscan(fid, '%d64 "%3s" %f %f %f %f*[\n]', 'delimiter', ',');
%     
% end
% 
% fclose(fid);

fileString = fixLineEndings(fileread(csv_file));
fileString = regexprep(fileString, '"[^"]*"','""');
Data = textscan(fileString,'%[^\n]',1,'HeaderLines',0,...
        'Delimiter',',','EndOfLine','\n');
headerLine = Data{1}{1};
hdr = strsplit(headerLine,',');
    
T = textscan(fileString,'%d64%3s%f%f%f%f','HeaderLines',1,...
        'Delimiter',',','EndOfLine','\n');

function str = fixLineEndings(str)
% Remove any \r\n, or \r with \n
% char(13) = '\r' and char(10) = '\n'

% make \r\n into \n
str(strfind(str,char([13 10]))) = '';
% make remaining \r into \n
str(str==char(13)) = char(10);

end

end

