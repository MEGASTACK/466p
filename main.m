function [ output_args ] = main( dat_path )
%MAIN Main script for segmentation pipeline
%   Run this to do the thing!
%   Matlab is alright.


%Load datas
[~, max_press] = pedo_extract(dat_path);


addpath JustinsGraphicalModelsToolboxPublic

% HOW TO: http://users.cecs.anu.edu.au/~jdomke/JGMT/
%
% features: cell array of structs (ie the samples)
% efeats: cell array of "edge features", use edgify_im
% loss_spec: 


feats = {};
efeats = [];
labels = [];
models = [];
loss_spec = [];
crf_type = 'linear_linear';

rho = 0;

options.viz         = @viz;
options.print_times = 0; % since this is so slow, print stuff to screen
options.gradual     = 1; % use gradual fitting
options.maxiter     = 1000;
options.rho         = rho;
options.reg         = 1e-4;
options.opt_display = 0;

%Do the training
params = train_crf(feats,efeats,labels,models,loss_spec,crf_type,options);

%Then test and evaluate


end

