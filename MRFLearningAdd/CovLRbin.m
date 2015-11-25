function C = CovLRbin(q,W)
% C = CovLRbin(q,W)
%
% input: q = p(s=1).
% W weight-matrix. Make sure there are zeros in the slots where there are no connections.
%
% output: C. The covariance matrix for biases and weights.
%
% AUTHOR: Max Welling 
 
 Nn = length(q);
 
 tmpl = (W~=0);
 z = sum(tmpl);
 [I,J] = find(triu(tmpl)==1);
 Ne = length(I);
  
 xi = Fxi2(W,q).*tmpl;
 
 dbdq = zeros(Nn);
 for i=1:Nn
     neighb_i = find(tmpl(i,:)==1);
     for k=1:Nn
         if k==i
             t1 = (1-z(i))/(q(i)*(1-q(i)));
             t2 = 0;
             for j=neighb_i
                 t2 = t2 + (1/(q(i)-xi(i,j)) + 1/(xi(i,j)+1-q(i)-q(j)));
             end
             dbdq(i,k) = t1 + t2;
         end
         if ismember(k,neighb_i)
             dbdq(i,k) = 1/(xi(i,k)+1-q(i)-q(k));
         end
     end
 end
 
 dbdxi = zeros(Nn,Ne);
 for i=1:Nn
     neighb_i = find(tmpl(i,:)==1);
     for b=1:Ne
         if (i==I(b) & ismember(J(b),neighb_i)) 
            dbdxi(i,b) = -(1/(q(i)-xi(I(b),J(b))) + 1/(xi(I(b),J(b))+1-q(i)-q(J(b))));
        end
        if (i==J(b) & ismember(I(b),neighb_i))
            dbdxi(i,b) = -(1/(q(i)-xi(J(b),I(b))) + 1/(xi(J(b),I(b))+1-q(i)-q(I(b))));
        end
    end
end

dWdq = zeros(Ne,Nn);
for a=1:Ne
    for k=1:Nn
        if I(a)==k 
            dWdq(a,k) = -(1/(xi(I(a),J(a))+1-q(I(a))-q(J(a))) + 1/(q(I(a))-xi(I(a),J(a))));
        end
        if J(a)==k
            dWdq(a,k) = -(1/(xi(I(a),J(a))+1-q(I(a))-q(J(a))) + 1/(q(J(a))-xi(I(a),J(a))));
        end
    end
end

dWdxi = zeros(Ne);
for a=1:Ne
    for b=1:Ne
        if I(a)==I(b) & J(a)==J(b)
            dWdxi(a,b) = 1/xi(I(a),J(a)) + 1/(xi(I(a),J(a))+1-q(I(a))-q(J(a))) + 1/(q(I(a))-xi(I(a),J(a)))...
                       + 1/(q(J(a))-xi(I(a),J(a)));
        end
    end
end

iC = [dbdq, dbdxi ; dWdq, dWdxi];

C = inv(iC);
        
        
        
        
        
        
        
          
         
         
     
             
 
 
