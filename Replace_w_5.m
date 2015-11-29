function [ new_matrix ] = Replace_w_5( old_matrix )
%Replace_w_5 replaces all of the non-zero values greater than 4 with a 5
%   Detailed explanation goes here
[k,t]=size(old_matrix);
for i = 1:k
    for j = 1:t
        if old_matrix(i,j)>4
            old_matrix(i,j)=5;
        end
    end
end
new_matrix = old_matrix;

end

