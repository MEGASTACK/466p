function [ Output ] = RegisterImage( Reference, Target )
%% RegisterImage
% Registers (aligns) the input image (Target) against the given Reference
% image. Requires/uses the MATLAB Image Processing Toolbox.

% Extend the image boundaries to avoid unwanted clipping.
[rows, cols] = size(Reference);
Reference = [zeros(rows, cols), Reference, zeros(rows, cols)];
Reference = [zeros(fix(rows/2), cols * 3); Reference; zeros(fix(rows/2), cols * 3)];
[rows, cols] = size(Target);
Target = [zeros(rows, cols), Target, zeros(rows, cols)];
Target = [zeros(fix(rows/2), cols * 3); Target; zeros(fix(rows/2), cols * 3)];

% Convert the input matricies to grayscale image formats.
% First, determine the maximum pressure value between the two images.
RefMaxVal = max(max(Reference));
TarMaxVal = max(max(Target));

% Convert both into grayscale images.
RefIm = mat2gray(Reference, [0, RefMaxVal]);
TarIm = mat2gray(Target,    [0, TarMaxVal]);

%imshowpair(RefIm, TarIm);
%pause;

% Perform registration
[optimizer, metric] = imregconfig('monomodal');
% Notes from Mathworks:
% Your registration results can improve if you adjust the optimizer or 
% metric settings. For example, if you increase the number of iterations in
% the optimizer, reduce the optimizer step size, or change the number of 
% samples in a stochastic metric, the registration improves to a point, at 
% the expense of performance.
optimizer.MaximumIterations = 300;
Output = imregister(TarIm, RefIm, 'rigid', optimizer, metric);

%imshowpair(RefIm, Output);

% The output is normalized to [0, 1], so first we can re-scale it back to
% its original format.
Output = Output .* TarMaxVal;

end

