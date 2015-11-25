function postSamples = HMCwrapper(data,A,priorHyperParams,numSamples,saveFile,leapFrogStep,varargin)

% parts of code segments (the actual HMC part) adapted from Max Welling's code
% Works only for 0/1 rep.
% Works only for grids.
%
% Add check for correct size etc. of varargin (to provide starting chain
% values)
%
% Add check if A infact represents a grid. 
%
% AUTHOR: Sridevi Parise

% check that data in 0/1 representation
if(~isempty( find(data~=0 & data~=1) ))
    fprintf(2,'Invalid Data. Accepts only 0/1 data\n');
    return;
end;

BURN_IN    = 0.5;
MAX_ITER   = 6000;
NUM_CHAINS = 5;
CONV_CHECK_ITER = 200;       % check convergence after every CONV_CHECK_ITER iter
CONV_SAVE_DIAG = 0;          % save MCMC diagnostics for possible visual plotting?
CONV_EPSILON = 0.1;
INIT_PRIOR_STD = 1;          % std. deviation for inital chain values (overdispersed?)
MAX_LAG        = 49;        % used to estimate thinning ( THIN_ITER can be atmost (1 + 2.MAX_LAG) )
ACORR_EPS      = 0.05;
PARTIAL_RES_INP = 1;
INIT_VAL_INP = 2;
NUM_SAMPLE_SETS = 1;         % while thinning save multiple sample sets if possible
N_LEAP_FROG   = 20;
LEAP_FROG_STEP = leapFrogStep;

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
adjModel.A = A;
adjModel.N = V;

rand('state',sum(100*(clock)));
randn('state',sum(100*(clock)));
%rand('state',3);
%randn('state',3);

% initialize chain
currModel.b(:,:) = zeros(V,NUM_CHAINS);
if( nargin==6 )              % no partial results or initial values
    maxSamplesStored = max( ceil( (1-BURN_IN)*MAX_ITER + BURN_IN*CONV_CHECK_ITER ), numSamples*(1+2*MAX_LAG)/NUM_CHAINS );
	samples = zeros(numParams,NUM_CHAINS,maxSamplesStored);
   	samples(:,:,1) = randn(numParams,NUM_CHAINS)*INIT_PRIOR_STD;
   	fprintf(1,'Metropolis: check for overdispersed initial values..\n');
	iter = 1;
	numRetSamples = 1;
	acceptNum = zeros(NUM_CHAINS,1);
    thinIter = 1;
    minIter = (numSamples)/(NUM_CHAINS);     % minimum number of iter required to get numSamples (assuming no thinning and no burn_in, updated later)
    converged = 0;
    convergedIter = MAX_ITER;
    maxLagAcorr = 1;
elseif( (nargin==8) & (varargin{1}==PARTIAL_RES_INP) )         % partial results provided
	parRes = load(varargin{2});
	samples = parRes.samples;
	iter = parRes.iter;
	numRetSamples = parRes.numRetSamples;	 
	acceptNum = parRes.acceptNum;
    thinIter = parRes.thinIter;
    minIter = parRes.minIter;
    converged = parRes.converged;
    convergedIter = parRes.convergedIter;
    maxLagAcorr = parRes.maxLagAcorr;
	clear parRes;
elseif( (nargin==8) & (varargin{1}==INIT_VAL_INP) )         % initial values provided
    maxSamplesStored = max( ceil( (1-BURN_IN)*MAX_ITER + BURN_IN*CONV_CHECK_ITER ), numSamples*(1+2*MAX_LAG)/NUM_CHAINS );
	samples = zeros(numParams,NUM_CHAINS,maxSamplesStored);
	samples(:,:,1) = varargin{2};
	iter = 1;
	numRetSamples = 1;
	acceptNum = zeros(NUM_CHAINS,1);
    thinIter = 1;
    minIter = (numSamples)/(NUM_CHAINS);     % minimum number of iter required to get numSamples (assuming no thinning and no burn_in, updated later)
    converged = 0;
    convergedIter = MAX_ITER;
    maxLagAcorr = 1;
else
	fprintf(1,'INCORRECT USAGE\n');
	return;
end;	
currModel.b(:,:) = samples(1:V,:,numRetSamples);
currModel.w = zeros(V,V,NUM_CHAINS);
for chain=1:NUM_CHAINS
    currModel.w(nzTriu_wIndMul(:,chain)) = samples(V+1:end,chain,numRetSamples);
    currModel.w(:,:,chain) = currModel.w(:,:,chain) + currModel.w(:,:,chain)';
end;

% compute log prob. and derivative of log prob. at current parameters
logOldP = zeros(NUM_CHAINS,1);
logOld_dP = zeros(numParams,NUM_CHAINS);
for chain = 1:NUM_CHAINS
    
    % model in adj. format
    adjModel.w = squeeze( currModel.w(:,:,chain) );
    adjModel.b = currModel.b(:,chain);
    % model in grid format
    gridModel = mapModelStructs( adjModel );
    
    logOldP(chain) = logPrior(priorHyperParams,samples(:,chain,numRetSamples));
    if( logOldP(chain)==1 )
        fprintf(2,'Error in logPrior\n');
        return;
    end;
    temp = (0.5*sum(sum(currModel.w(:,:,chain).*sum_xi_xj))) + (sum(currModel.b(:,chain).*sum_xi));
    logOldP(chain) = logOldP(chain) + temp;
    sNodes = runJTgrid( gridModel );
    logZ   = findLogZ( sNodes );       % compute log(Z) using junction tree marginals
    logOldP(chain) = logOldP(chain) - (N * logZ);
    
    % the derivative
    marginals = findJTmarginals( sNodes );
    logOld_dP(:,chain) = [sum_xi;sum_xi_xj(nzTriu_wInd)] - (N*[marginals.node;marginals.edge(nzTriu_wInd)]) - (priorHyperParams.invCov*(samples(:,chain,numRetSamples)-priorHyperParams.mu));
    
end;

logNewP = logOldP;
logNew_dP = logOld_dP;
while(1)
    
    if( mod(iter,500)==0 )
        fprintf(1,'HMC ITER = %d\n',iter);
        save(saveFile,'samples','numRetSamples','iter','acceptNum','thinIter','minIter','converged','convergedIter','maxLagAcorr');
    end;
    
    % run one step of HMC 
    potParams = samples(:,:,numRetSamples);
    hmc_p = randn(numParams,NUM_CHAINS);   % initialize momenta
    K_x = sum(hmc_p.^2)/2;  % compute kinetic energy
    U_x = -logOldP';
    Dx = -logOld_dP;
    
    for i=1:N_LEAP_FROG   % loop over leapfrog steps
        hmc_p = hmc_p - 0.5*LEAP_FROG_STEP * Dx;
        potParams = potParams + LEAP_FROG_STEP * hmc_p;       
        % find derivatives        
        potModel.b = potParams(1:V,:);
        potModel.w = zeros(V,V,NUM_CHAINS);
        for chain = 1:NUM_CHAINS
            potModel.w(nzTriu_wIndMul(:,chain)) = potParams(V+1:end,chain);
            potModel.w(:,:,chain) = potModel.w(:,:,chain) + potModel.w(:,:,chain)';
        end;
        for chain = 1:NUM_CHAINS
            % model in adj. format
            adjModel.w = squeeze( potModel.w(:,:,chain) );
            adjModel.b = potModel.b(:,chain);
            % model in grid format
            gridModel = mapModelStructs( adjModel );
            sNodes = runJTgrid( gridModel );
            marginals = findJTmarginals( sNodes );
            logNew_dP(:,chain) = [sum_xi;sum_xi_xj(nzTriu_wInd)] - (N*[marginals.node;marginals.edge(nzTriu_wInd)]) - (priorHyperParams.invCov*(potParams(:,chain)-priorHyperParams.mu));
        end;
        % end find derivatives
        Dt =  -logNew_dP;
        hmc_p = hmc_p - 0.5*LEAP_FROG_STEP*Dt;
        Dx = Dt;
    end
    
    % find logP at proposal        
    potModel.b = potParams(1:V,:);
    potModel.w = zeros(V,V,NUM_CHAINS);
    for chain = 1:NUM_CHAINS
        potModel.w(nzTriu_wIndMul(:,chain)) = potParams(V+1:end,chain);
        potModel.w(:,:,chain) = potModel.w(:,:,chain) + potModel.w(:,:,chain)';
    end;
    for chain = 1:NUM_CHAINS
        % model in adj. format
        adjModel.w = squeeze( potModel.w(:,:,chain) );
        adjModel.b = potModel.b(:,chain);
        % model in grid format
        gridModel = mapModelStructs( adjModel );
        logNewP(chain) = logPrior(priorHyperParams,potParams(:,chain));
        if( logNewP(chain)==1 )
            fprintf(2,'Error in logPrior\n');
            return;
        end;
        temp = (0.5*sum(sum(potModel.w(:,:,chain).*sum_xi_xj))) + (sum(potModel.b(:,chain).*sum_xi));
        logNewP(chain) = logNewP(chain) + temp;
        sNodes = runJTgrid( gridModel );
        logZ   = findLogZ( sNodes );       % compute log(Z) using junction tree marginals
        logNewP(chain) = logNewP(chain) - (N * logZ);
    end;
    % end find logP
    U_t = -logNewP'; % use your own energy function here
    K_t = sum(hmc_p.^2)/2;
    
    DH = U_x + K_x - U_t - K_t; % compute difference total energy
    
    hmc_P = min(1,exp(DH));
    accept = rand(1,NUM_CHAINS) < hmc_P;  % accept according to metroplis-hastings criterium
    acceptVec(iter,:) = accept;
    samples(:,:,numRetSamples+1) = samples(:,:,numRetSamples);
    samples(:,accept,numRetSamples+1) = potParams(:,accept);
    
    logOldP(accept) = logNewP(accept);
    logOld_dP(:,accept) = logNew_dP(:,accept);
    acceptNum = acceptNum + accept';
    iter = iter+1;
    numRetSamples = numRetSamples+1; 
    % end run one step of HMC
    
    if( mod(iter,CONV_CHECK_ITER)==0 )
        
        % test convergence
        if( ~converged )
            % remove additional burn-in samples at this check pt.
            temp = floor( BURN_IN*CONV_CHECK_ITER );        % samples to remove (valid only at conv. chekc pt.s)
            samples(:,:,1:numRetSamples-temp) = samples(:,:,temp+1:numRetSamples);
            numRetSamples = numRetSamples - temp;
            % end remove

            % compute R_c
			tempSamp = shiftdim( samples(:,:,1:numRetSamples), 2 );
            [Rmpsrf,neff_mpsrf,Vmpsrf,Wmpsrf,Bmpsrf] = mpsrf( tempSamp );
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
                %minIter = (minIter*thinIter) + floor(convergedIter*BURN_IN);          % update minIter
                minIter = (minIter*max([thinIter NUM_SAMPLE_SETS])) + floor(convergedIter*BURN_IN);
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
    %minIter = (minIter*thinIter) + floor(convergedIter*BURN_IN); 
    minIter = (minIter*max([thinIter NUM_SAMPLE_SETS])) + floor(convergedIter*BURN_IN);
end;

% collect samples from all chains
%sampleSetsPossible = min(thinIter,NUM_SAMPLE_SETS);
sampleSetsPossible = NUM_SAMPLE_SETS;
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
    fprintf(1,'Metropolis: MAX_ITER reached. samples less than requested\n');
end;
if( ~converged )
    fprintf(1,'Metropolis: MAX_ITER reached. chains not converged \n');
end;
if( maxLagAcorr>ACORR_EPS )
    fprintf(1,'Metropolis: MAX_LAG reached at max acorr=%0.2f. possible insufficient thinning\n',maxLagAcorr);
end;


acceptRate = acceptNum/iter;
fprintf(1,'Metropolis: total iter=%d mean acceptRate=%0.2f thinIter=%d\n',iter,mean(acceptRate),thinIter);
save acceptVec.mat acceptVec -V6;
if( CONV_SAVE_DIAG )
    save convDiagV6.mat convDiag acceptRate;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes a gaussian prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function logP = logPrior(hyperParams,params)

if( isempty(hyperParams) )
    logP = 0;
    return;
end;

% check compatibility of sizes
numParams = length(params);
[r,c] = size(hyperParams.cov);
if( (length(hyperParams.mu)~=numParams) | (r~=c) | (r~=numParams) )
    fprintf(2,'Error in gaussian prior - size mismatch\n');
    logP = 1;
    return;
end;

logP = -0.5*(params-hyperParams.mu)'*(hyperParams.invCov)*(params-hyperParams.mu);
logP = logP - (0.5*numParams*log(2*pi));
logP = logP - (0.5*log(hyperParams.detCov));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes log(Z) using the fact that: Z = 1/p(0) for a boltzmann m/c
% p(0) is computed using the marginals got by running JT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function logZ = findLogZ( sNodes )

num_sNodes  = size(sNodes.nodes,1);
endBy2plus1 = size(sNodes.pot,2)/2 + 1;     % 2^R + 1  


log_p0 = log( sNodes.pot(1,1) );     % contribution from first supernode/clique variables of JT
% contribution from remaining variables (one from each clique)
for clq=2:num_sNodes
    log_p0 = log_p0 + log( sNodes.pot(clq,1)/( sNodes.pot(clq,1)+sNodes.pot(clq,endBy2plus1) ) );
end;

logZ = -log_p0;




