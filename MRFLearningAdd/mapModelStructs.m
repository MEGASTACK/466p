function outModel = mapModelStructs( inModel ) 

% Converts between the two model struct formats:
% 1. The adjacency format with fields:
%    (assumes node numbering is same as column-wise on the corressponding grid)
%    (assume valid adjacency matrix)
%
% 2. The grid format with fields:
%
% AUTHOR: Sridevi Parise

if( isfield(inModel,'A') )  % input format is ADJ
    R = 1;
    while( (R<inModel.N) & (inModel.A(R,R+1)) ) % find numRows
        R = R+1;
    end;
    C = (inModel.N)/R;
    
    outModel.numRows = R;
    outModel.numCols = C;
    outModel.alpha = reshape(inModel.b,R,C);
    for r=1:R
        for c=1:C-1
            outModel.wHor(r,c) = inModel.w(R*c-R+r,R*c+r);
        end;
    end;
    for r=1:R-1
        for c=1:C
            outModel.wVer(r,c) = inModel.w(R*c-R+r,R*c-R+r+1);
        end;
    end;
    
    if( isfield(inModel,'Vind') )   % if visible node info present
        outModel.Vind = inModel.Vind;
        outModel.Hind = inModel.Hind;
    end;
else       % input format is GRID
    R = inModel.numRows;
    C = inModel.numCols;
    
    outModel.N = R*C;
    outModel.b = reshape(inModel.alpha,outModel.N,1);
    outModel.A = zeros(outModel.N,outModel.N);
    outModel.w = zeros(outModel.N,outModel.N);
    for r=1:R
        for c=1:C-1
            outModel.w(R*c-R+r,R*c+r) = inModel.wHor(r,c);
            outModel.A(R*c-R+r,R*c+r) = 1;
        end;
    end;
    for r=1:R-1
        for c=1:C
            outModel.w(R*c-R+r,R*c-R+r+1) = inModel.wVer(r,c);
            outModel.A(R*c-R+r,R*c-R+r+1) = 1;
        end;
    end;
    outModel.A = outModel.A + outModel.A';
    outModel.w = outModel.w + outModel.w';
    
    if( isfield(inModel,'Vind') )   % if visible node info present
        outModel.Vind = inModel.Vind;
        outModel.Hind = inModel.Hind;
    end;
end;


