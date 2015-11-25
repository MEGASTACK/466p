function F = Fbethe2(W,ta,xi,q)
% function F = Fbethe2(W,ta,xi,q)
% computes the bethe free energy.
% Author: Max Welling
 

 [D,N] = size(q);
 
 tmpl = (W~=0); 
 z = sum(tmpl,2);
 TMPL = repmat(tmpl,[1,1,N]);
 Z = repmat(z,[1,N]);
 
 I = repmat(eye(D),[1,1,N]);;
 qq1 = permute(repmat(q,[1,1,D]),[1,3,2]);
 qq2 = permute(qq1,[2,1,3]);

 ns = 1e-100;
 W = repmat(W,[1,1,N]);
 
 
 
 F = - squeeze(sum(sum((xi.*W).*TMPL)))'/2 - sum(ta.*q) + ...
       squeeze(sum(sum(((xi).*log(xi+I+ns)).*TMPL)))'/2 +...
       squeeze(sum(sum(((xi+1-qq1-qq2).*log(xi+1-qq1-qq2+ns)).*TMPL)))'/2+...
       squeeze(sum(sum(((qq1-xi).*log(qq1-xi+ns)).*TMPL)))'/2 +...
       squeeze(sum(sum(((qq2-xi).*log(qq2-xi+ns)).*TMPL)))'/2 +...
       sum((1-Z).*(q.*log(q+ns)+(1-q).*log(1-q+ns)));
      
  
 F=real(F);

 
