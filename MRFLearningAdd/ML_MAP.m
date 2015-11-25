function [model] = ML_MAP(data,R,C,priorCov)

% Estimate the MAP parameters for a rectangular grid boltzmann m/c,
% using the junction tree algorithm for exact inference
% Inputs:
%           data : (#samples X #nodes) matrix containing training samples
%                  Each node can take values from {+1,-1} or {0,1}
%			R    : # rows in the grid
%           C    : # columns in the grid
%        priorCov: (numParams X numParams) covariance matrix defining a zero
%                mean gaussian prior over the parameters. While indexing parameters 
%                biases apppear first, then wHor and then wVer (all ordered column wise). 
%                [] if no prior (ML estimation)
%
% Returns:
%           model: struct with fields
%                  numRows: R
%                  numCols: C
%                  alpha  : (R X C) matrix of node biases (MAP estimates)
%                  wHor   : (R X C-1) matrix of horizontal edge weights (MAP estimates)
%                  wVer   : (R-1 X C) matrix of vertical edge weights (MAP estimates)
%
% The node value representation ( +1/-1 or 0/1) intended by the user is guessed from the training samples.       
% If using -1/+1, the data is first mapped to 0/1 and the model params are learned. 
% These learned params are then mapped back to the -1/1 case.
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

% if -1/1 rep. then map data to 0/1
if( ~REP_01 )
    data = 0.5*( data+1 );
end;
    
	
% now learn params

MAX_ITERATIONS = 50000;
EPSILON        = 0.00001;
EPSILON1       = 0.0001;
rho            = 10;   % the initial step size
MONITOR        = 0;
rand('state',sum(100*(clock)));
randn('state',sum(100*(clock)));
%rand('state',3);
%randn('state',3);

N = size(data,1);
if( ~isempty(priorCov) )
    USE_PRIOR = 1;
    priorCovInv = inv(priorCov);
    priorGrad_alpha = zeros(R,C);
    priorGrad_wHor = zeros(R,C-1);
    priorGrad_wVer = zeros(R-1,C);
else
    USE_PRIOR = 0;
end;    

% extract the sufficient statistics from data
sStat.alpha = reshape((1/N)*sum(data,1),[R C]);
for r=1:R
    for c=1:C-1
        sStat.wHor(r,c) = (1/N)*sum(data(:,sub2ind([R C],r,c)).*data(:,sub2ind([R C],r,c+1)));
    end;
end;
for r=1:R-1
    for c=1:C
        sStat.wVer(r,c) = (1/N)*sum(data(:,sub2ind([R C],r,c)).*data(:,sub2ind([R C],r+1,c)));
    end;
end;
% end extract sufficient statistics

%initialize the model
model.numRows = R;
model.numCols = C;
model.alpha = randn(R,C);
model.wHor  = randn(R,C-1);
model.wVer  = randn(R-1,C);

if(MONITOR)  % save initial params
    paramVec.alpha = zeros(MAX_ITERATIONS,R,C);
    paramVec.wVer = zeros(MAX_ITERATIONS,R-1,C);
    paramVec.wHor = zeros(MAX_ITERATIONS,R,C-1);
    paramVec.alpha(1,:,:) = model.alpha;
    paramVec.wVer(1,:,:) = model.wVer;
    paramVec.wHor(1,:,:) = model.wHor;
end;

% run JT inference ( results used for computing likelihood and gradient )
sNodes = runJTgrid(model);

% compute prior contribution to gradient at current params
if(USE_PRIOR)
    priorGrad = -priorCovInv*[model.alpha(:);model.wHor(:);model.wVer(:)];
    priorGrad_alpha(:) = priorGrad(1:R*C);
    priorGrad_wHor(:) = priorGrad(R*C+1:R*(2*C-1));
    priorGrad_wVer(:) = priorGrad(R*(2*C-1)+1:end);
else
    priorGrad_alpha = 0;
    priorGrad_wHor = 0;
    priorGrad_wVer = 0;
end;

if(USE_PRIOR)
    lPos = logL(data,sNodes) - 0.5*([model.alpha(:);model.wHor(:);model.wVer(:)]'*priorCovInv*[model.alpha(:);model.wHor(:);model.wVer(:)]); % log posterior
else
    lPos = logL(data,sNodes);
end;
%llVec = lPos;
iter = 0;
newModel = model;

updates = 0;
while(iter<MAX_ITERATIONS)
    
    if( mod(iter,100)==0 )
        fprintf(1,'ITER=%d LL=%e rho=%f\n',iter,lPos,rho);
    end;
    iter = iter + 1;

    marginals = findGridMarginals( sNodes );       % find marginals required for gradient computation
    newModel.alpha = model.alpha + ( rho*( sStat.alpha-marginals.node+(priorGrad_alpha/N) ) );    % update alpha
    newModel.wHor = model.wHor + ( rho*( sStat.wHor-marginals.eHor+(priorGrad_wHor/N) ) );       % update wHor
    newModel.wVer = model.wVer + ( rho*( sStat.wVer-marginals.eVer+(priorGrad_wVer/N) ) );       % update wVer
     
    % run JT inference with new model params
    new_sNodes = runJTgrid(newModel);
    
    if(USE_PRIOR)
        newlPos = logL(data,new_sNodes)- 0.5*([newModel.alpha(:);newModel.wHor(:);newModel.wVer(:)]'*priorCovInv*[newModel.alpha(:);newModel.wHor(:);newModel.wVer(:)]);
    else
        newlPos = logL(data,new_sNodes);
    end;
    
    if( newlPos<lPos )
        rho = rho/2;
        continue;
    end;
    
    updates = updates + 1;
    rho = (1.1)*rho;
    
    if(MONITOR) % save params at this iteration
        paramVec.alpha(updates+1,:,:)    = newModel.alpha;
        paramVec.wVer(updates+1,:,:)    = newModel.wVer;
        paramVec.wHor(updates+1,:,:)    = newModel.wHor;
    end;
    
    if( (max(abs(newModel.alpha(:)-model.alpha(:)))<EPSILON1) & (max(abs(newModel.wVer(:)-model.wVer(:)))<EPSILON1) & (max(abs(newModel.wHor(:)-model.wHor(:)))<EPSILON1) )
        model = newModel;
        sNodes = new_sNodes;
        lPos = newlPos;
        break;
    end;
    
    model = newModel;
    sNodes = new_sNodes;
    
    % compute prior contribution to gradient at new params
    if(USE_PRIOR)
        priorGrad = -priorCovInv*[model.alpha(:);model.wHor(:);model.wVer(:)];
        priorGrad_alpha(:) = priorGrad(1:R*C);
        priorGrad_wHor(:) = priorGrad(R*C+1:R*(2*C-1));
        priorGrad_wVer(:) = priorGrad(R*(2*C-1)+1:end);
    else
        priorGrad_alpha = 0;
        priorGrad_wHor = 0;
        priorGrad_wVer = 0;
    end;
    
    
%     if( mod(iter,100)==0 )            % possible schedule to adapt rho
%         rho = rho/2;
%     end;
    
%     %if( (abs((newll-ll)/ll)<EPSILON) & (rho>MAX_RHO) )
%     if( abs((newll-ll)/ll)<EPSILON )
%         ll = newll;
%         llVec = [llVec ll];
%         break;
%     end;
    
    lPos = newlPos;
%   llVec = [llVec lPos];
end;


% find equivalent params in {-1/1} rep. if necessary
if( ~REP_01 )
    neighbors = getNeighbors(R,C);
    for i=1:(R*C)       % ( maybe more efficient ways to do this with a neater design )
        if( (~isempty(neighbors{i}.nodesHor))&(~isempty(neighbors{i}.nodesVer)) )
            model.alpha(i) = 0.5*model.alpha(i) + 0.25*sum(model.wHor(neighbors{i}.eHor)) + 0.25*sum(model.wVer(neighbors{i}.eVer));
        elseif( isempty(neighbors{i}.nodesHor) )
            model.alpha(i) = 0.5*model.alpha(i) + 0.25*sum(model.wVer(neighbors{i}.eVer));
        else
            model.alpha(i) = 0.5*model.alpha(i) + 0.25*sum(model.wHor(neighbors{i}.eHor));
        end;
    end;
    model.wHor = 0.25*model.wHor;
    model.wVer = 0.25*model.wVer;
end;
% end map params

fprintf(1,'Iterations used = %d MAX_ITER=%d updates=%d\n',iter,MAX_ITERATIONS,updates);

if( MONITOR )
    save fooV6.mat paramVec llVec -V6;
end;

return;    % MLgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ll = logL(data,sNodes)

N = size(data,1);
num_sNodes = size(sNodes.nodes,1);
R = size(sNodes.nodes,2)-1;
C = (num_sNodes/R) + 1;
powVec = cumprod([1 2*ones(1,R)])';    % [1 2 2^2 2^3 ....]

ll = 0;

% first sNode (first R+1 variables)
ind = 1 + data(:,1:R+1)*powVec;        % get linear indices into the potential function for all data samples
ll = ll + sum( log(sNodes.pot(1,ind)) ); % contributions from all samples for these variables

% remaining sNodes (remaining variables, one in each sNode)
for clq=2:num_sNodes
    
    ind = 1 + [data(:,clq:clq+R-1) zeros(N,1)]*powVec;
    p0 = sNodes.pot(clq,ind);
    ind = 1 + [data(:,clq:clq+R-1) ones(N,1)]*powVec;
    p1 = sNodes.pot(clq,ind);
    % normalize
    norm = p0+p1;
    p0 = p0./norm;
    p1 = p1./norm;
    
    ll = ll + log(p0)*( 1-data(:,clq+R) ) + log(p1)*data(:,clq+R);  % add p0 if variable takes value 0 etc. 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
