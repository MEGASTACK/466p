function [ score ] = JaccardScore( A, B )
%% JaccardScore
% Takes in 2 matricies and computes the Jaccard Score of the matricies A
% and B. Strictly speaking, it calculates the Jaccard Score based on
% element-by-element matching.

Intersection = A == B;
Intersection = sum(sum(Intersection));

score = Intersection / (numel(A) + numel(B) - Intersection);

end

