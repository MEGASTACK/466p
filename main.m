function [ output_args ] = main( dat_path )
%MAIN Main script for segmentation pipeline
%   Run this to do the thing!
%   Matlab is alright.


%Load datas
[~, max_press] = pedo_extract(dat_path);


efeats = [];
labels = [];
models = [];
loss_spec = [];
crf_type = 'linear_linear';

N     = 1;  % number of training images
siz   = 50; % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 2;  % this problem is binary (switch me)..


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

model_hash = repmat({[]},1000,1000);
fprintf('building models...\n')
for n=1:N
    [ly lx lz] = size(ims{n});
    if isempty(model_hash{ly,lx});
        model_hash{ly,lx} = gridmodel(ly,lx,nvals);
    end
end
models = cell(N,1);
for n=1:N
    [ly lx lz] = size(ims{n});
    models{n} = model_hash{ly,lx};
end



edge_params = {{'const'},{'diffthresh'},{'pairtypes'}};
fprintf('computing edge features...\n')
efeats = cell(N,1);
parfor n=1:N
    efeats{n} = edgeify_im(ims{n},edge_params,models{n}.pairs,models{n}.pairtype);
end



fprintf('splitting data into a training and a test set...\n')
k = 1;
% [who_train who_test] = kfold_sets(N,5,k);
% 
% ims_train     = ims(who_train);
% feats_train   = feats(who_train);
% efeats_train  = efeats(who_train);
% labels_train  = labels(who_train);
% labels0_train = labels0(who_train);
% models_train  = models(who_train);
% 
% ims_test     = ims(who_test);
% feats_test   = feats(who_test);
% efeats_test  = efeats(who_test);
% labels_test  = labels(who_test);
% labels0_test = labels0(who_test);
% models_test  = models(who_test);

ims_train = ims;
ims_test = ims;
feats_train = feats;
feats_test = feats;
efeats_train = efeats;
efeats_test = efeats;
labels_train = labels;
labels_test = labels;
labels0_train = labels0;
labels0_test = labels0;
models_train = models;
models_test = models;


    % visualization function
    function viz(b_i)
        % here, b_i is a cell array of size nvals x nvars
        M = 5;
        for n=1:M
            [ly lx lz] = size(ims_train{n});
            subplot(3,M,n    ); miximshow(reshape(b_i{n}',ly,lx,nvals),nvals);
            subplot(3,M,n+  M); imshow(ims_train{n})
            subplot(3,M,n+2*M); miximshow(reshape(labels_train{n},ly,lx),nvals);

        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end

loss_spec = 'trunc_cl_trwpll_5';


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

