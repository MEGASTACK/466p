function model = PL_MAP(data,A,priorCov)

% Estimate the PL parameters (MAP sense) for a boltzmann m/c
% Inputs: 
%          data: (N X V) matrix containing the data samples. 
%                 N is the total number of samples and V is the total number of nodes                                         
%                 Each node can take values from {+1,-1} or {0,1}.
%             A: (V X V) adjacency matrix defining the graph structure
%      priorCov: (numParams X numParams) covariance matrix defining a zero
%                mean gaussian prior over the parameters. While indexing parameters 
%                biases apppear first then weights (ordered column wise in the upper triangular part
%                of A, no main diagonal). 
%                [] if no prior (ML sense estimation)
% Returns:
%          model: (1X1) struct array with fields
%                   N: the number of nodes
%                   A: adjacency matrix                 
%                   b: the biases  ( PL estimates )
%                   w: the edge weights ( PL estimates )
%
% The node value representation ( +1/-1 or 0/1) intended by the user is detected from the training samples.       
% If using 0/1, the data is first mapped to +1/-1 and the model params are learned. 
% These learned params are then mapped back to the 0/1 case.
%
%
% AUTHOR: Sridevi Parise


% find node value representation ( assuming error free data )
minVal = min(min(data));
if( minVal==-1 )
    REP_01 = 0;
elseif( minVal==0 )
    REP_01 = 1;
else
    fprintf(2,'Ambiguous-all ones dataset\n');
    return;
end;
  
% if 0/1 rep. then map data to +1/-1
if( REP_01 )
    data = 2*data - 1;
end;
    
% now learn params

MAX_ITERATIONS = 10000;
ADAPT_PTS        = [500:500:10000];
RHO_CHANGE       = 0.5*ones(1,length(ADAPT_PTS));
EPSILON = 0.00001;
EPSILON1 = 0.0001;
momentum = 0.9;
rho            = 0.5;   % the initial step size
batchSize = 400;
MONITOR = 1;
rand('state',sum(100*(clock)));
randn('state',sum(100*(clock)));

N = size(data,1);
V = size(data,2);
numAdaptPts = length( ADAPT_PTS );
if( ~isempty(priorCov) )
    USE_PRIOR = 1;
    priorCovInv = inv(priorCov);
else
    USE_PRIOR = 0;
end;    
nzTriu_wInd = find( triu(A,1) );       % linear indices for non-zero wt.s (or edges) above main diagonal

if( N<batchSize )
    batchInd = [0 N];  % indices of batch boundaries
else
    batchInd = [0:batchSize:N];
    if( batchInd(end)<N )
        batchInd = [batchInd N];
    end;
end;
numBatches = length(batchInd)-1;

% extract some batch statistics from data
for i=1:numBatches
    batchData = data(batchInd(i)+1:batchInd(i+1),:);
    batchLen  = size(batchData,1);
    sStat(i).sum_xi    = sum(batchData,1);
    sStat(i).sum_xi_xj = batchData'*batchData;
    sStat(i).sum_xi_mat = repmat(sStat(i).sum_xi,V,1);
end;
% end extract

%initialize the model
model.N = V;
model.A = A;
model.b = randn(V,1);
model.w = randn(V,V).*A;
model.w = 0.5 * (model.w + model.w');    % make it symmetric


if(MONITOR)  % save initial params
    nonZero_wInd = find(model.w);
    paramVec.b = zeros(MAX_ITERATIONS+1,V);
    paramVec.w = zeros(MAX_ITERATIONS+1,length(nonZero_wInd));
    paramVec.b(1,:) = model.b;
    paramVec.w(1,:) = model.w(nonZero_wInd);
end;

%initialize the delta's (for adding momentum)
delta_b = zeros(V,1);
delta_w = zeros(V,V);

%initialize local rho's (adaptive local learning rates)
rho_b = ones(V,1);
rho_w = ones(V,V);

iter = 0;
adapt = 1;
while(iter<MAX_ITERATIONS)
    
    iter = iter + 1;
    if( mod(iter,100)==0 )
        fprintf(1,'ITER=%d\n',iter);
    end;
    
    currBatch = mod(iter-1,numBatches)+1;
    batchData = data(batchInd(currBatch)+1:batchInd(currBatch+1),:);
    batchLen  = size(batchData,1);
    
    prevModel = model;   % save current model
    
    % compute prior contribution to gradient
    if(USE_PRIOR)
        priorGrad = -priorCovInv*[model.b;model.w(nzTriu_wInd)];
        priorGrad_b = priorGrad(1:V);
        priorGrad_w = zeros(V,V);
        priorGrad_w(nzTriu_wInd) = priorGrad(V+1:end);
        priorGrad_w = priorGrad_w + priorGrad_w';        % make it symmetric
    else
        priorGrad_b = 0;
        priorGrad_w = 0;
    end;

    % compute \sum_{n_i}w_{i,n_i}x_{n_i}^(r) where n_i are neighbors of x_i for each i and data point r
    condStat = batchData*model.w;
    
    temp = (2)./( 1+exp( -2*(repmat(model.b,1,batchLen)+condStat') ) );
    % update b's
    delta_b = ((rho/batchLen)*rho_b).*( sStat(currBatch).sum_xi' + batchLen - sum(temp,2) + priorGrad_b ) + (momentum*delta_b);
    model.b = model.b + delta_b;
    
    % update w's
    temp1 = temp*batchData;
    delta_w = ((rho/batchLen)*rho_w).*( 2*sStat(currBatch).sum_xi_xj+sStat(currBatch).sum_xi_mat+sStat(currBatch).sum_xi_mat'-temp1-temp1'+priorGrad_w ).*model.A + (momentum*delta_w);
    model.w = model.w + delta_w;
    
    if(MONITOR) % save params at this iteration
        paramVec.b(iter+1,:)    = model.b;
        paramVec.w(iter+1,:)    = model.w(nonZero_wInd);
    end;
    
%     if( mod(iter,200)==0 )            % possible schedule to adapt rho
%         rho = rho/2;
%     end;
    
    if( adapt <= numAdaptPts)
        if( iter==ADAPT_PTS(adapt) )
            rho = rho*RHO_CHANGE(adapt);
            adapt = adapt+1;
        end;
    end;
    
    if( (max(abs(model.b-prevModel.b))<EPSILON1) & (max(abs(model.w(:)-prevModel.w(:)))<EPSILON1) )
        break;
    end;
    
end;

% map params if necessary
if( REP_01 )
    model.b = 2*model.b - 2*sum(model.w,2);
    model.w = 4*model.w;
end;   

fprintf(1,'Iterations used = %d MAX_ITER=%d\n',iter,MAX_ITERATIONS);

if(MONITOR) % visualize the changes (in +1/-1 rep.)
    plotParamVar( paramVec );
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = plotParamVar( paramVec )

save fooPLv6.mat paramVec -V6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
