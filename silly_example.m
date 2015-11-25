function [ ] = silly_example( )
%MAIN Main script for segmentation pipeline
%   Run this to do the thing!
%   Matlab is alright.






% First, set up the parameters of the problem
N     = 5;  % number of training images
siz   = 50; % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 2;  % this problem is binary
% Next, generate the graph for this CRF. (A simple pairwise grid)
model = gridmodel(siz,siz,nvals);
% Now, generate the data. Basically, we make noisy images, then smooth them to make the true (discrete) output values, and then add noise to make the input.

for n=1:N
    x{n} = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric'));
    t = rand(size(x{n}));
    noiselevel = 1.25;
    y{n} = x{n}.*(1-t.^noiselevel) + (1-x{n}).*t.^noiselevel;
end
% In this case, the data looks like this. First the true labels x, and then the noisy inputs y.

% Next, we make the features and the labels. The features consist of simply a constant of one, and the input image y itself.

for n=1:N
    feats{n}  = [y{n}(:) 1+0*x{n}(:)];
    labels{n} = x{n}+1;
end

% no edge features here
efeats = []; % none
% Next, a visualization function. This takes a cell array of predicted beliefs as input, and shows them to the screen during training. This is totally optional, but very useful if you want to understand what is happening in your training run.

    % visualization function
    function viz(b_i)
        % here, b_i is a cell array of size nvals X nvars
        for n=1:N
            subplot(3,N,n    ); imshow(reshape(b_i{n}(2,:),siz,siz));
            subplot(3,N,n+  N); imshow(reshape(feats{n}(:,1),siz,siz));
            subplot(3,N,n+2*N); imshow(reshape(labels{n}-1,siz,siz));

        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end
% Next, we pick a string to specify the loss and inference method. In this case, we choose truncated fitting with the clique logistic loss based on TRW with five iterations.

loss_spec = 'trunc_cl_trw_5';
% Other options include 'pert_ul_trw_1e5' (perturbation, univariate logistic loss, TRW, threshold of 1e-5) 'em_mnf_1e5' (Surrogate Expectation-Maximization based on mean-field with a threshold of 1e-5 (simplifies to surrogate likelihood with no hidden variables) or 'trunc_em_trwpll_10' (Truncated surrogate EM based on multithreaded TRW with 10 iterations).
% Next, we set up some parameters for the training optimization.

crf_type  = 'linear_linear';
options.derivative_check = 'off';
options.viz         = @viz;
options.rho         = rho;
options.print_times = 1;
options.nvals       = nvals;
% Finally we actually optimize.
p = train_crf(feats,efeats,labels,model,loss_spec,crf_type,options)


%Then test and evaluate


end

