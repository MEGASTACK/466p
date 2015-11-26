function [ Output ] = RegisterImage( Reference, Target )
%% RegisterImage
% Registers (aligns) the input image (Target) agains the given Reference
% image. Requires/uses the MATLAB Image Processing Toolbox.

% Convert the input matricies to grayscale image formats.
% First, determine the maximum pressure value between the two images.
MaxValue = max(max([Target, Reference]));

% Convert both into grayscale images.
RefIm = mat2gray(Reference, [0, MaxValue]);
TarIm = mat2gray(Target,    [0, MaxValue]);

% Perform registration
[optimizer, metric] = imregconfig('multimodal');
Output = imregister(TarIm, RefIm, 'rigid', optimizer, metric);

% The output is normalized to [0, 1], so first we can re-scale it back to
% its original format.
Output = Output .* MaxValue;

end

