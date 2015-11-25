function postSamples = briefLangevin(data,A,priorHyperParams,numSamples,K,saveFile,IS_GRID,varargin)

% Works only for 0/1 rep.
%
% Add check for correct size etc. of varargin (to provide starting chain
% values)
%
%
% AUTHOR: Sridevi Parise

% check that data in 0/1 representation
if(~isempty( find(data~=0 & data~=1) ))
    fprintf(2,'Invalid Data. Accepts only 0/1 data\n');
    return;
end;

BURN_IN = 0.5;           % fraction of chain to discard
STEP_SZ = 0.01;          % step-size for langevin method 
MAX_ITER   = 400000;
NUM_CHAINS = 5;
CONV_CHECK_ITER = 200;       % check convergence after every CONV_CHECK_ITER iter
CONV_SAVE_DIAG = 0;          % save MCMC diagnostics for possible visual plotting?
CONV_EPSILON = 0.1;
INIT_PRIOR_STD = 1;          % std. deviation for inital chain values (overdispersed?)
MAX_LAG        = 49;        % used to estimate thinning ( THIN_ITER can be atmost (1 + 2.MAX_LAG) )
ACORR_EPS      = 0.05;
PARTIAL_RES_INP = 1;
INIT_VAL_INP = 2;
NUM_SAMPLE_SETS = 20;         % while thinning save multiple sample sets if possible


N = size(data,1);
V = size(data,2);
nzTriu_wInd = find( triu(A) );      % indicies of non-zero weights (each edge counted only once) (into 2D array)
nzTriu_wIndMul = [];                % indicies of non-zero weights for case of multiple chains (into 3D array)
for chain = 1:NUM_CHAINS
    nzTriu_wIndMul(:,chain) = nzTriu_wInd+(chain-1)*V*V;
end;
numParams = V + length(nzTriu_wInd);
sum_xi_xj = data'*data;
sum_xi = sum(data)';

if( IS_GRID )
    fprintf(1,'Assuming A represents a valid grid (with node numbering column-wise on corresponding grid)\n');
    % find a 2-coloring of the grid
    dim1 = 1;
    while( A(dim1,dim1+1)==1 )
        dim1 = dim1 + 1;
    end;
    dim2 = size(A,1)/dim1;
    color1Ind = [];
    for i=1:dim2
        startInd = (i-1)*dim1 + 1;
        endInd = startInd+dim1-1;
        if( mod(i,2)==0 )
            color1Ind = [color1Ind startInd:2:endInd];
        else
            color1Ind = [color1Ind startInd+1:2:endInd];
        end;
    end;
    color2Ind = setdiff( [1:size(A,1)], color1Ind );
    % end find 2-coloring
    color1Ind = color1Ind';
    color2Ind = color2Ind';
    clear dim1 dim2 startInd endInd i;
end;

rand('state',sum(100*(clock)));
randn('state',sum(100*(clock)));
%rand('state',3);
%randn('state',3);

% initialize chains
currModel.b(:,:) = zeros(V,NUM_CHAINS);
if( nargin==7 )            % no partial results or intial values
    maxSamplesStored = max( ceil( (1-BURN_IN)*MAX_ITER + BURN_IN*CONV_CHECK_ITER ), numSamples*(1+2*MAX_LAG)/NUM_CHAINS );
	samples = zeros(numParams,NUM_CHAINS,maxSamplesStored);
    samples(:,:,1) = randn(numParams,NUM_CHAINS)*INIT_PRIOR_STD;
    fprintf(1,'briefLangevin: check for overdispersed initial values..\n');
	iter = 1;
	numRetSamples = 1;
    thinIter = 1;
    minIter = (numSamples)/(NUM_CHAINS);     % minimum number of iter required to get numSamples (assuming no thinning and no burn_in, updated later)
    converged = 0;
    convergedIter = MAX_ITER;
    maxLagAcorr = 1;
elseif( (nargin==9) & (varargin{1}==PARTIAL_RES_INP) )    % partial results provided
	parRes = load(varargin{2});
	samples = parRes.samples;
	iter = parRes.iter;
	numRetSamples = parRes.numRetSamples;	
    thinIter = parRes.thinIter;
    minIter = parRes.minIter;
    converged = parRes.converged;
    convergedIter = parRes.convergedIter;
    maxLagAcorr = parRes.maxLagAcorr;
	clear parRes;
elseif( (nargin==9) & (varargin{1}==INIT_VAL_INP) )        % initial values provided
    maxSamplesStored = max( ceil( (1-BURN_IN)*MAX_ITER + BURN_IN*CONV_CHECK_ITER ), numSamples*(1+2*MAX_LAG)/NUM_CHAINS );
	samples = zeros(numParams,NUM_CHAINS,maxSamplesStored);
    samples(:,:,1) = varargin{2};
	iter = 1;
	numRetSamples = 1;
    thinIter = 1;
    minIter = (numSamples)/(NUM_CHAINS);     % minimum number of iter required to get numSamples (assuming no thinning and no burn_in, updated later)
    converged = 0;
    convergedIter = MAX_ITER;
    maxLagAcorr = 1;
else
	fprintf(1,'briefLangevin: INCORRECT USAGE\n');
    return;
end;
	
while(1)
    
    if( mod(iter,5000)==0 )
        fprintf(1,'Langevin ITER = %d\n',iter);
        save(saveFile,'samples','numRetSamples','iter','thinIter','minIter','converged','convergedIter','maxLagAcorr');
    end;
    
    % convert to struct format
    currModel.N = V;
    currModel.b(:,:) = samples(1:V,:,numRetSamples);
    currModel.w = zeros(V,V,NUM_CHAINS);
    for chain=1:NUM_CHAINS
        currModel.w(nzTriu_wIndMul(:,chain)) = samples(V+1:end,chain,numRetSamples);
        currModel.w(:,:,chain) = currModel.w(:,:,chain) + currModel.w(:,:,chain)';
    end;
   
    
    % find K-step reconstructions using current model
    kStepSamples = zeros(N,V,NUM_CHAINS);
    for chain=1:NUM_CHAINS
        currChainModel.N = currModel.N;
        currChainModel.b = currModel.b(:,chain);
        currChainModel.w = currModel.w(:,:,chain);
        if( IS_GRID )
            kStepSamples(:,:,chain) = KstepMCMCgrid(data,currChainModel,K,color1Ind,color2Ind);
        else
            kStepSamples(:,:,chain) = KstepMCMC(data,currChainModel,K);
        end;
    end;
    
    
    % likelihood contribution to gradient
    for chain = 1:NUM_CHAINS
        llGrad_w(:,:,chain) = sum_xi_xj - (kStepSamples(:,:,chain)'*kStepSamples(:,:,chain));  
        llGrad_b(:,chain) = sum_xi - sum(kStepSamples(:,:,chain))';
        llGrad(:,chain) = [llGrad_b(:,chain); llGrad_w(nzTriu_wIndMul(:,chain))];
    end;    
    
    % prior contribution to gradient
    priorGrad = findPriorGrad(priorHyperParams,samples(:,:,numRetSamples));
    
    samples(:,:,numRetSamples+1) = samples(:,:,numRetSamples) + ( (STEP_SZ^2)*(0.5)*(llGrad+priorGrad) ) + ( STEP_SZ*randn(numParams,NUM_CHAINS) );
    
    iter = iter+1;
    numRetSamples = numRetSamples+1;
    
    if( mod(iter,CONV_CHECK_ITER)==0 )
        
        % test convergence
        if( ~converged )
        	% remove additional burn-in samples at this check pt.
        	temp = floor( (BURN_IN)*CONV_CHECK_ITER );        % samples to remove (valid only at conv. chekc pt.s)
        	samples(:,:,1:numRetSamples-temp) = samples(:,:,temp+1:numRetSamples);
        	numRetSamples = numRetSamples - temp;
        	% end remove

            % compute R_c
            [Rmpsrf,neff_mpsrf,Vmpsrf,Wmpsrf,Bmpsrf] = mpsrf( shiftdim( samples(:,:,1:numRetSamples), 2 ) );
            fprintf(1,'R=%f, W=%f, B=%f\n',Rmpsrf,rcond(Wmpsrf),rcond(Bmpsrf));
            if( CONV_SAVE_DIAG )
                convDiag.R(iter/CONV_CHECK_ITER) = Rmpsrf;
                convDiag.V(iter/CONV_CHECK_ITER) = det(Vmpsrf);
                convDiag.W(iter/CONV_CHECK_ITER) = det(Wmpsrf);
            end;
            % end compute R_c
          
            % check if R_c close to 1 and estimate thinning if converged
            if( Rmpsrf<(1+CONV_EPSILON) )
                converged = 1;
                convergedIter = iter;
                % estimate thinning iterations required
                acorrTimes = zeros(NUM_CHAINS,numParams);
                for chain = 1:NUM_CHAINS
                    acorrVec = acorr( squeeze(samples( :,chain,1:numRetSamples ))', MAX_LAG );
                    acorrTimes(chain,:) = 1 + 2*sum(acorrVec);
                    maxLagAcorr(chain) = max( acorrVec(end,:) );      % acorr at MAX_LAG (max over all params)
                end;
                thinIter = ceil( max( mean( acorrTimes ) ) );
                maxLagAcorr = max( maxLagAcorr );          % max over all chains
                % end estimate thinning
                minIter = (minIter*thinIter) + floor(convergedIter*BURN_IN);          % update minIter
            end;
        end;
        
    end;
    % end check pt.
    
    if( (iter==MAX_ITER) | ((iter>=minIter)&converged) )
        break;
    end;
    
end;
% end MCMC

if( ~converged )
    %estimate thinning even though not converged
    acorrTimes = zeros(NUM_CHAINS,numParams);
    for chain = 1:NUM_CHAINS
        acorrVec = acorr( squeeze(samples(:,chain,1:numRetSamples))', MAX_LAG );
        acorrTimes(chain,:) = 1 + 2*sum(acorrVec);
        maxLagAcorr(chain) = max( acorrVec(end,:) );      % acorr at MAX_LAG (max over all params)
    end;
    thinIter = ceil( max( mean( acorrTimes ) ) );
    maxLagAcorr = max( maxLagAcorr );          % max over all chains
    % end estimate thinning
    minIter = (minIter*thinIter) + floor(convergedIter*BURN_IN); 
end;

% collect samples from all chains
sampleSetsPossible = min(thinIter,NUM_SAMPLE_SETS);
postSamples = [];
for s=1:sampleSetsPossible
    postSamples(s).samples = [];
end;
totIter = iter;
if( ~converged )
    startIter = ( numRetSamples - ceil((1-BURN_IN)*iter) ) + 1;
else
    startIter = 1;
end; 
for s=1:sampleSetsPossible
    for chain = 1:NUM_CHAINS
        postSamples(s).samples = [postSamples(s).samples squeeze( samples(:,chain,startIter+s-1:thinIter:numRetSamples) )];
    end;
    if( size(postSamples(s).samples,2)>numSamples )
        postSamples(s).samples = postSamples(s).samples(:,end-numSamples+1:end);
    end;
end; 
% end collect samples

% print possible failures
if( iter<minIter )
    fprintf(1,'briefLangevin: MAX_ITER reached. samples less than requested\n');
end;
if( ~converged )
    fprintf(1,'briefLangevin: MAX_ITER reached. chains not converged \n');
end;
if( maxLagAcorr>ACORR_EPS )
    fprintf(1,'briefLangevin: MAX_LAG reached at max acorr=%0.2f. possible insufficient thinning\n',maxLagAcorr);
end;

fprintf(1,'briefLangevin: total iter=%d thinIter=%d\n',iter,thinIter);
if( CONV_SAVE_DIAG )
    save convDiagV6.mat convDiag -V6;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes gradient for a gaussian prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function priorGrad = findPriorGrad(hyperParams,params)

% check compatibility of sizes
numParams = size(params,1);
[r,c] = size(hyperParams.cov);
if( (length(hyperParams.mu)~=numParams) | (r~=c) | (r~=numParams) )
    fprintf(2,'Error in gaussian prior - size mismatch,[..Ignoring prior..]\n');
    priorGrad = 0;
    return;
end;

priorGrad = -(hyperParams.invCov*params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samples = KstepMCMC(data,model,K)

N = size(data,1);

samples = data;
for k=1:K
    for i=1:model.N       % sample model.N nodes randomly, given the rest
        node = ceil( rand(1,1)*model.N );   % node to sample
        % find p(x_node=1|x_{-node}) for each sample
        p = 1./( 1 + exp( -(model.b(node)+samples*sparse(model.w(node,:)')) ) );
        % sample node
        r = rand(N,1);
        samples(:,node) = (r<p);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samples = KstepMCMCgrid(data,model,K,color1Ind,color2Ind)

N = size(data,1);
bMat = repmat( model.b',N,1 );
numColor1 = length(color1Ind);
numColor2 = length(color2Ind);

samples = data;
for k=1:K
    
    % sample color1 nodes, given the rest
    p = 1./( 1 + exp( -(bMat(:,color1Ind)+samples*sparse(model.w(:,color1Ind))) ) );
    r = rand(N,numColor1);
    samples(:,color1Ind) = (r<p);
    
    % sample color2 nodes, given the rest
    p = 1./( 1 + exp( -(bMat(:,color2Ind)+samples*sparse(model.w(:,color2Ind))) ) );
    r = rand(N,numColor2);
    samples(:,color2Ind) = (r<p);
    
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



