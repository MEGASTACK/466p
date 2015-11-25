function [posMu,posCov] = posDistMLgrid(D,R,C,priorCov)
%
%
%
% AUTHOR: Sridevi Parise

N = size(D,1);
V = size(D,2);

% find approximate MAP model
fprintf(1,'Finding MAP value.....\n');
modelMAP = ML_MAP(D,R,C,priorCov);
modelMAP = mapModelStructs( modelMAP );          % convert to adj. format
nzTriu_wInd = find( triu(modelMAP.A) );
lambdaMAP = [modelMAP.b;modelMAP.w(nzTriu_wInd)];   % the parameter vector

% find marginals for MAP model using BP
fprintf(1,'Finding marginals using BP.....\n');
[p_xi1,p_xi1_xj1,BPiter] = fastBPbin(modelMAP.w,modelMAP.b);
fprintf(1,'BPiter = %d\n',BPiter);

% find feature covariance using BP-LR
C = CovLRbin(p_xi1,modelMAP.w);

% compute the posterior approximation
% sum_xi = sum(D,1)';
% sum_xi_xj = D'*D;
% sum_fi = [sum_xi;sum_xi_xj(nzTriu_wInd)];  % x_i and x_i_x_j are features
% E_fi = [p_xi1;p_xi1_xj1(nzTriu_wInd)];

% modelMAPgrid = mapModelStructs( modelMAP );
% sNodes = runJTgrid( modelMAPgrid );
% marginals = findJTmarginals( sNodes );
% E_fi = [marginals.node;marginals.edge(nzTriu_wInd)];

invPriorCov = inv(priorCov);
posCov = (1/N) * inv( C + (invPriorCov/N) );
% posMu = posCov*( sum_fi-(N*E_fi)-(invPriorCov*lambdaMAP) );
% posMu = posMu + lambdaMAP;
posMu = lambdaMAP;
