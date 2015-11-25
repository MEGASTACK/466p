function [posMu,posCov] = posDist(D,A,priorCov)
%
%
%
% AUTHOR: Sridevi Parise

N = size(D,1);
V = size(D,2);
nzTriu_wInd = find( triu(A) );

% find approximate MAP model
fprintf(1,'Finding MAP value.....\n');
%modelMAP = PL_MAP(D,A,priorCov);
if( isValidData(D,A) )          % no empty cells
    initModel = PMM(D,A);       % initialize with PMM
    modelMAP = BP_MAP(D,A,priorCov,initModel);
else
    modelMAP = BP_MAP(D,A,priorCov);
end;
lambdaMAP = [modelMAP.b;modelMAP.w(nzTriu_wInd)];   % the parameter vector

% find marginals for MAP model using BP
fprintf(1,'Finding marginals using BP.....\n');
[p_xi1,p_xi1_xj1,BPiter] = fastBPbin(modelMAP.w,modelMAP.b);
fprintf(1,'BPiter = %d\n',BPiter);

% find feature covariance using BP-LR
C = CovLRbin(p_xi1,modelMAP.w);

% compute the posterior approximation
sum_xi = sum(D,1)';
sum_xi_xj = D'*D;
sum_fi = [sum_xi;sum_xi_xj(nzTriu_wInd)];  % x_i and x_i_x_j are features
E_fi = [p_xi1;p_xi1_xj1(nzTriu_wInd)];

invPriorCov = inv(priorCov);
posCov = (1/N) * inv( C + (invPriorCov/N) );
%posMu = posCov*( sum_fi-(N*E_fi)-(invPriorCov*lambdaMAP) );
%posMu = posMu + lambdaMAP;
posMu = lambdaMAP;
