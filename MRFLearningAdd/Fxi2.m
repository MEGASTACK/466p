function xi = Fxi2(W,q);
%  function xi = Fxi2(W,q);
% input: weight-matrix W and marginals q --> p(s=1)
% output: xi is correlations between neighbors: p(s_i=1,ps_j=1).
%
% AUTHOR: Max Welling
 
 N = length(q);
 
 Qi = repmat(q',N,1);
 Qj = repmat(q, 1,N);
 
 tmpl = (W ~= 0);

 
 bt = 1./(exp(W + ~tmpl) - 1);
 
 qiqj = q*q';
 QQ = (bt + Qi + Qj);
 
 xi1 = qiqj .* ~tmpl;
  
 xi2 = ((QQ - sign(bt).*sqrt(QQ.^2 - 4*(1+bt).*qiqj)) / 2) .* tmpl;
 
 xi = real((xi1 + xi2) .* ~eye(N));

 
