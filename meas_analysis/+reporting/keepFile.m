function keep = testFile(meta, freq, location_str, loc_logic)
% Filter out this file
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 4
    loc_logic = true;
end

if nargin < 3
    location_str = NaN;
end

if nargin < 2
    error 'frequency specification is required'
end

if nargin <1
    error 'meta file not provided'
end

keep = false;

% filter
if ~isnan(freq)
    if meta.Frequency_GHz_num == freq
        keep = true;
    else
        keep = false;
        return;
    end
end

if ~isnan(location_str)
    is_loc = ~isempty(strfind(meta.MatFile_str,location_str));
    if (is_loc && loc_logic) || (~is_loc && ~loc_logic)
        keep = true;
    else
        keep = false;
        return;
    end
end

end