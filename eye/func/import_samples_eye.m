function [ T, hdr ] = import_samples_eye( csv_file )

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

