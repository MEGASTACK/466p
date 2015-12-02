function [ images ] = RegularizeData( filenames )
%% RegularizeData
%   Input a cell of size n by 1 of filenames (lst pedograph files).
%   Output an x by y by n matrix of the max pressure maps, where each page
%   in the 3rd dimension is one of the input files, padded so that it's the
%   same size in x and y as the largest image.

%% Extract the max pressure maps from all of the files
rawMaxMaps = cell(length(filenames), 1);
maxSize = [0, 0];
for i = 1 : length(filenames)
    [~, maxPressure, rows, cols] = pedo_extract(filenames{i});
    rawMaxMaps{i} = maxPressure;
    maxSize = max(maxSize, [rows, cols]);
end

%% Pad each image and add it to the output.
images = zeros([maxSize, length(rawMaxMaps)]);
for i = 1 : length(rawMaxMaps)
    
    padIm = rawMaxMaps{i};
    [m, n] = size(padIm);
    
    if any(lt([m, n], maxSize))
        dM = maxSize(1) - m;
        dN = maxSize(2) - n;
        
        % Add dM rows, dN cols
        padIm = [padIm, zeros(m, dN)];
        padIm = [padIm; zeros(dM, n+dN)];
    end
    
    images(:,:,i) = padIm;
end

end

