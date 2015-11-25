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


feat_params = {{'patches',0},{'position',1},{'fourier',1},{'hog',8}};


fprintf('loading data and computing feature maps...\n');
N = 1;
for n=1:N
    % load data
%     lab = importdata([labdir lab_names(n).name]);
    lab = max_press;
    im  = double(max_press)/255;
    ims{n}  = im;
    labels0{n} = max(0,lab+1);

    % compute features
    feats{n}  = featurize_im(ims{n},feat_params);

    % reduce resolution for speed
%     ims{n}    = imresize(ims{n}   ,rez,'bilinear');
%     feats{n}  = imresize(feats{n} ,rez,'bilinear');
%     labels{n} = imresize(labels0{n},rez,'nearest');

    % reshape features
    [ly lx lz] = size(feats{n});
    feats{n} = reshape(feats{n},ly*lx,lz);
end



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

