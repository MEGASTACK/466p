function [ index ] = RandIndex( X, Y )
%% RandIndex
% Computes the rand index between the matricies A and B. Assumes that X and
% Y are 2 different labelings of S.

% Flatten X and Y
X = X(:);
Y = Y(:);

% https://en.wikipedia.org/wiki/Rand_index
% Given a set of n elements S and two partitions of S to compare, X: a
% partition of S into r subsets and Y: a partition of S into s subsets and:
% a: The number of pairs of elements in S that are in the same set in X and
%    the same set in Y.
a = 0;
% b: the number of pairs of elements in S that are in different sets in X
%    and different sets in Y.
b = 0;

% Compute a and b
for i = 1 : numel(X)
    for j = 1 : numel(X)
        if i == j
            continue;
        end
        
        if X(i) == Y(j)
            a = a + 1;
        else
            b = b + 1;
        end
    end
end

% The Rand index is:
% R = a + b / n choose 2
index = (a + b) / nchoosek(numel(X), 2);

end

