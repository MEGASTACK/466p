function [ data ] = LoadPedoData(rootDataDir, normalize)
%% LoadPedoData
%   Loads all max-maps from lst files in the given directory as an NxMxP
%   matrix, where each page P in the matrix is one of the max maps.
%       If rootDataDir is not specified, then it defaults to 'FinalData'
%       If the input <normalize> is exactly 1, then the data is also 
%       normalized to the range [0, 1].

% Apply defaults
if nargin < 2
    normalize = 0;
end
if nargin < 1
    rootDataDir = 'FinalData';
end
    
% Collect filenames
fNames = GetLstFilenames(rootDataDir);

% Get the data from the list of filenames
data = RegularizeData(fNames);

% Normalize if we need to (translate into grayscale image)
if normalize == 1
    for i = 1 : size(data, 3)
        raw = data(:,:,i);
        data(:,:,i) = mat2gray(raw, [0, max(max(raw))]);
    end
end

end

