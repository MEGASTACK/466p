function [q,xi,loop] = fastBPbin(W,b)
% function [q,xi,loop] = fastBPbin(W,b)
% input: W: weights-matrix and b: biases.
% output: q marginals p(s=1)
%         xi correlations between neighbors p(s_i=1,ps_j=1).
%
%  AUTHOR: Max Welling
%
%  5/31/05 ---- possible adhoc modifications - sridevi


MAX_LOOP = 10000;

[D,N] = size(b);
if N>1
   B = squeeze(permute(repmat(b,[1,1,D]),[1,3,2]));
elseif N==1
   B = repmat(b,[1,D]);
end

  tmpl = (W~=0);
  Tmpl = squeeze(repmat(tmpl,[1,1,N]));

  Mess = squeeze(zeros(D,D,N).*Tmpl);    
  
  al = squeeze(repmat((exp(W)-1),[1,1,N]));
    
  threshold = 1e-10;
  change = threshold+1;
  Nloop = 1000;
  loop = 1;
 
  DAMP = [0*ones(1,100) 0.1*ones(1,100) 0.2*ones(1,100) 0.3*ones(1,100) 0.4*ones(1,100) ...
        0.5*ones(1,100) 0.6*ones(1,100) 0.7*ones(1,100) 0.8*ones(1,100) 0.9*ones(1,100)];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 while ( abs(change)>threshold & loop<=MAX_LOOP ) 
     
   if( loop<Nloop )   
       damp = DAMP(loop); 
   else
       damp = DAMP(Nloop);
   end;
   
   if N > 1
      SMess = squeeze(sum(Mess,1));
      SMess = permute(repmat(SMess,[1,1,D]),[3,1,2]);
      SMess = permute((SMess - Mess),[2,1,3]);
      
      Sigm = 1./(1+exp( -B - SMess ));
      
      NMess = ( damp * Mess + (1-damp) * log( 1 + al.*Sigm ) ) .* Tmpl ;
      
   elseif N==1
      SMess = sum(Mess,1);
      SMess = repmat(SMess,[D,1]);
      SMess = (SMess - Mess)';
      
      NMess = damp * Mess + (1-damp) * (log( 1 + al.*sigmoid(B + SMess) ) .* tmpl);
      
   end
   
   change = max(max(max(abs(NMess - Mess))));
   
  % fprintf('ITERATION: %.3d, CHANGE: %.3u \n',loop, change )
  
   Mess = NMess;
   loop = loop + 1;
   
   end %while %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
      
   if N > 1
      q = exp(squeeze(sum(Mess,1)) + b);
      q = q ./ (1 + q);
      
      xi11 = exp( repmat(W,[1,1,N]) + B + permute(B,[2,1,3]) + SMess + permute(SMess,[2,1,3]) );
      xi12 = exp( B + SMess );
      xi21 = exp( permute(B,[2,1,3]) + permute(SMess,[2,1,3]) );
      Sxi = 1 + xi11 + xi12 + xi21;
      xi = (xi11 ./ Sxi) .* Tmpl;
      
   elseif N == 1
      q = exp(sum(Mess,1)' + b);
      q = q ./ (1 + q);

      xi11 = exp( W + B + B' + SMess + SMess');
      xi12 = exp( B + SMess );
      xi21 = exp( B' + SMess' );
      Sxi = 1 + xi11 + xi12 + xi21;
      xi = (xi11 ./ Sxi) .* tmpl;
   end
    
      
     
   
   
   
   
   
   
   
   
      
    
      
