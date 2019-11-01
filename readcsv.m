function [ D, hdr ] = readcsv( filename )
% How hard is this

fid = fopen(filename, 'r');
D = textscan(fid,'%s%d','Delimiter',',');
fclose(fid);

if nargout > 1
   hdr = D(1,:);
   D = D(2:end,:);
end


end

