% l2_norm
%
% Computes the l2_norm of a vector. If given a N x M numeric array, this
% will return a N x 1 vector of the l2_norms of each row.
%
% Implementation by Jesse Hellemn at Rice University, Fall 2016
function norm = l2_norm(V)
    norm = sqrt(sum(V .^ 2, 2));
end
