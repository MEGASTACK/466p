function model = BP_MAP(data,A,priorCov,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters (MAP sense) for a boltzmann m/c using
% Loopy BP Learning
% INPUTS: 
%        data        : (N X V) matrix where N = # data points and V = # visible nodes
%					    Each variable can take values from {+1,-1} or {0,1}				
%        A           : Adjacency matrix defining the graph structure for
%                      the model to learn 
%        priorCov    : (numParams X numParams) covariance matrix defining a zero
%                       mean gaussian prior over the parameters. While indexing parameters 
%                       biases apppear first then weights (ordered column wise in the upper triangular part
%                       of A, no main diagonal). 
%                       [] if no prior (ML sense estimation)
%        varargin     : can provide optional initial model 
% RETURNS:
%          model: (1X1) struct array with fields
%                   N: the number of nodes
%                   A: adjacency matrix                 
%                   b: the biases  ( BP estimates )
%                   w: the edge weights ( BP estimates )
%
% The node value representation ( +1/-1 or 0/1) intended by the user is guessed from the training samples.       
% If using +1/-1, the data is first mapped to 0/1 and the model params are learned. 
% These learned params are then mapped back to the +1/-1 case.
%
%
% AUTHOR : Sridevi Parise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
% if +1/-1 rep. then map data to 0/1
if( ~REP_01 )
    data = 0.5*( data+1 );
end;
    
% now learn params

MAX_ITERATIONS   = 2000;
ADAPT_PTS        = [2000];
RHO_CHANGE       = 0.5*ones(1,length(ADAPT_PTS));
EPSILON1         = 0.00001;
rho              = 0.5;   % the initial step size
ADAPT_INT        = 500;     % adapt rates after every ADAPT_INT iterations
SMOOTH_WIN       = 50;       % take mean over a window of size SMOOTH_WIN 
LOCAL_ADAPT      = 0; 
momentum         = 0.9;
MONITOR          = 0;     % monitor the change in params? (useful for testing if above learning parameter settings are reasonable)
rand('state',sum(100*(clock)));
randn('state',sum(100*(clock)));
% rand('state',110);
% randn('state',2000);

numAdaptPts = length( ADAPT_PTS );
if( ~isempty(priorCov) )
    USE_PRIOR = 1;
    priorCovInv = inv(priorCov);
else
    USE_PRIOR = 0;
end;    
nzTriu_wInd = find( triu(A,1) );       % linear indices for non-zero wt.s (or edges) above main diagonal


N = size(data,1);
V = size(data,2);      % number of nodes

% extract statistics from data
sStat.E_xi    = sum(data,1)'/N;
sStat.E_xi_xj = (data'*data)/N;
% end extract

%initialize the model
if( nargin==4 )
    model = varargin{1};
elseif( nargin==3 )
    model.N = V;
    model.A = A;
    model.b = randn(model.N,1)*0.001;
    model.w = (randn(model.N,model.N)*0.001).*A;
    model.w = 0.5 * (model.w + model.w');    % make it symmetric
else
    fprintf(1,'BP_MAP: Incorrect Usage\n');
end;

nonZero_wInd = find(model.w); 
numEdges = length(nonZero_wInd);
if(MONITOR)  % save initial params
    paramVec.b = zeros(MAX_ITERATIONS+1,model.N);
    paramVec.w = zeros(MAX_ITERATIONS+1,length(nonZero_wInd));
    paramVec.b(1,:)    = model.b;
    if( numEdges~=0 )
        paramVec.w(1,:) = model.w(nonZero_wInd);
    end;
end;

%initialize the delta's (for adding momentum)
delta_b = zeros(model.N,1);
delta_w = zeros(model.N,model.N);

%initialize local rho's (adaptive local learning rates)
rho_b = ones(model.N,1);
rho_w = ones(model.N,model.N);

iter = 0;
bpIter = 0;
adapt = 1;
while(iter<MAX_ITERATIONS)
    
    if(mod(iter,100)==0)
        fprintf(1,'ITER=%d BP_ITER=%d\n',iter,bpIter);
    end;
    iter = iter + 1;
    
	prevModel = model;             % save current model
    
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


    % expected values from equilibrium distribution using LBP
    [E_xi,E_xi_xj,bpIter] = fastBPbin(model.w,model.b);
    
     % update b
     delta_b = (rho*( rho_b.*(sStat.E_xi-E_xi+(priorGrad_b/N)) )) + (momentum*delta_b);  % the current delta
     model.b = model.b + delta_b;
         
    % update w
    delta_w =  (rho*( rho_w.*((sStat.E_xi_xj-E_xi_xj+(priorGrad_w/N)).*model.A) )) + (momentum*delta_w);
    model.w = model.w + delta_w;
        
    if(MONITOR) % save params at this iteration
        paramVec.b(iter+1,:)    = model.b;
        if( numEdges~=0 )
            paramVec.w(iter+1,:)  = model.w(nonZero_wInd);
        end;
    end;
    
    % possible local learning rate adaption scheme
    if( LOCAL_ADAPT & mod(iter,ADAPT_INT)==0 )
       
        mean1.b = mean( paramVec.b(iter-ADAPT_INT+1:iter-ADAPT_INT+SMOOTH_WIN,:) ); 
        mean1.w = mean( paramVec.w(iter-ADAPT_INT+1:iter-ADAPT_INT+SMOOTH_WIN,:) );
        mean2.b = mean( paramVec.b(iter-(ADAPT_INT/2)+1:iter-(ADAPT_INT/2)+SMOOTH_WIN,:) ); 
        mean2.w = mean( paramVec.w(iter-(ADAPT_INT/2)+1:iter-(ADAPT_INT/2)+SMOOTH_WIN,:) );
        mean3.b = mean( paramVec.b(iter-SMOOTH_WIN+1:iter,:) ); 
        mean3.w = mean( paramVec.w(iter-SMOOTH_WIN+1:iter,:) );
        signPrev.b = sign(mean2.b-mean1.b);
        signPrev.w = sign(mean2.w-mean1.w);
        signCurr.b = sign(mean3.b-mean2.b);
        signCurr.w = sign(mean3.w-mean2.w);
        
        % if sign not same, decrease rate
        temp = 1 - ( 0.4*ones(model.N,1).*(signPrev.b.*signCurr.b==-1)' );
        rho_b = rho_b.*temp;
        temp = 1 - ( 0.4*ones(numEdges,1).*(signPrev.w.*signCurr.w==-1)' );
        rho_w(nonZero_wInd) = rho_w(nonZero_wInd).*temp;
        % if sign same, increase rate  
        temp = 1 + ( 0.2*ones(model.N,1).*(signPrev.b.*signCurr.b==1)' );
        rho_b = rho_b.*temp;
        temp = 1 + ( 0.2*ones(numEdges,1).*(signPrev.w.*signCurr.w==1)' );
        rho_w(nonZero_wInd) = rho_w(nonZero_wInd).*temp;
        
    end;
              
        
    %if( mod(iter,500)==0 )            % possible schedule to adapt rho
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
if( ~REP_01 )
    model.b = (0.5)*model.b + (0.25)*sum(model.w,2);
    model.w = (0.25)*model.w;
end;

if(MONITOR) % visualize the changes (in 0/1 rep.)
    plotParamVar( paramVec );
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = plotParamVar( paramVec )

save fooBP.mat paramVec -V6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
