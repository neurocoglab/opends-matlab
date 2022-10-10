function result = compare_versions(ver1, varargin)
% Compares two version numbers.
%
% Example: 7.4.0.287
%   major release.minor release.patch level.build
%   Everything except the major release is optional.
% 
% Arguments:
%     ver1 - Version number one
%     ver2 - Version number two; default: version cmd output
% 
% Return values:
%     result - Index of part which differs:
%              0 if equal
%            > 0 if ver1 > ver2
%            < 0 if ver1 < ver2
%            NaN if ver1 | ver2 are not a version number
% 
% Author: Johannes Kissel
% Last modified: 29.01.2008

    % check args
    if isempty(varargin) % no ver2 given
        [s, e, i, ver2] = matchver(version);
    else
        ver2 = varargin{1};
    end

    if ~isver(ver1) | ~isver(ver2) % ver1 | ver2 are not a version number
        result = NaN;
        return
    end

    % init result
    result = 0;

    % precalcs
    ver1 = sscanf(ver1, '%d.');
    ver2 = sscanf(ver2, '%d.');

    lV1 = length(ver1);
    lV2 = length(ver2);

    % loop over parts of ver1
    for i = 1:lV1
        if i > lV2 % ver1 longer than ver2
            result = i;
            return
        else
            nV1 = ver1(i);
            nV2 = ver2(i);

            if nV1 > nV2 % ver1 > ver2
                result =  i;
                return
            end

            if nV1 < nV2 % ver1 < ver2
                result = -i;
                return
            end
        end
    end

    if lV1 < lV2 % ver1 shorter than ver2
        result = -(lV1 + 1);
    end
end

function result = isver(in)
    result = 1;

    if ~ischar(in)
        result = 0;
    elseif length(strfind(in, '.')) > 3
        result = 0;
    else
        [s, e] = matchver(in);

        if isempty(s) | s > 1 | e < length(in)
            result = 0;
        end
    end
end

function varargout = matchver(in)
    [s, e, i, m] = regexp(in, '\d+(\.\d+)*', 'once');
    varargout = {s, e, i, m};
end

