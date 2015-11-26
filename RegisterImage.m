function [ Output ] = RegisterImage( Reference, Target )
%% RegisterImage
% Registers (aligns) the input image (Target) against the given Reference
% image. Requires/uses the MATLAB Image Processing Toolbox.

% Convert the input matricies to grayscale image formats.
% First, determine the maximum pressure value between the two images.
MaxValue = max(max([Target, Reference]));

% Convert both into grayscale images.
RefIm = mat2gray(Reference, [0, MaxValue]);
TarIm = mat2gray(Target,    [0, MaxValue]);

%imshowpair(RefIm, TarIm, 'diff');
%pause

% Perform registration
[optimizer, metric] = imregconfig('multimodal');
% Notes from Mathworks:
% Your registration results can improve if you adjust the optimizer or 
% metric settings. For example, if you increase the number of iterations in
% the optimizer, reduce the optimizer step size, or change the number of 
% samples in a stochastic metric, the registration improves to a point, at 
% the expense of performance.
%optimizer.MaximumIterations = 300;
Output = imregister(TarIm, RefIm, 'rigid', optimizer, metric);

%imshowpair(RefIm, Output, 'diff');

% The output is normalized to [0, 1], so first we can re-scale it back to
% its original format.
Output = Output .* MaxValue;

end

