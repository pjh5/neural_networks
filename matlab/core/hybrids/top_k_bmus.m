% top_k_bmus
%
% Given a N x D matrix of N data vectors and an M x D matrix of M prototype
% weight vectors, this returns an N x M matrix X_mapped where X_mapped[i, j] =
% k_weights[k] if prototype[k,:] is the kth BMU of X[i,:]. If the length of
% k_weights is k, then each row of the returned matrix will have exactly k
% non-zero elements at the indices corresponding to the k BMUs of that row.
%
% PARAMETERS
%
% prototypes    An M x D matrix of M prototype weight vectors, each of
%               dimension D, arranged into M rows.
%
% X             An N x D matrix of N data vectors, each of dimension D,
%               arranged in rows.
%
% k_weights     A vector of dimension K, these are the values that will be
%               inserted into the returned matrix. For example, if k_weights =
%               [0.6, 0.3, 0.1], then X_mapped[i, BMU(i)] will be 0.6,
%               X_mapped[i, 2nd_BMU(i)] will be 0.3, and X_mapped[i,
%               3rd_BMU(i)] will be 0.1. These k_weights should be the "weight"
%               that the kth BMU should have.
%
% RETURNS
%
% X_mapped      An N x M matrix. Rows correspond to data vectors, columns
%               correspond to prototypes. For each row i, the top K BMUs of
%               data vector i will be non-zero, and every other element will be
%               0. The kth BMU of data vector i will be the kth element of
%               k_weights.
%
% Written by Jesse Hellemn at Rice University, 2016
function [X_mapped] = top_k_bmus(prototypes, X, k_weights)

    % Cache number of BMUs to take
    n_k = numel(k_weights);

    % Allocate memory for mapped X matrix
    X_mapped = zeros(size(X, 1), size(prototypes, 1));

    % TODO find some variation of rowfun that works on matrices to do the
    % following:
    % rowfun(@(x) som_winner(protos, x), cell2table(mat2cell(X, ones(M,1))))
    % Until then loop
    for xi = 1:size(X, 1)

        % Calculate distances to every prototype
        data_dists = l2_norm(bsxfun(@minus, X(xi, :), prototypes));

        % Take K BMUs
        for ki = 1:n_k
            [~, w_idx] = min(data_dists);
            X_mapped(xi, w_idx) = k_weights(ki);
            data_dists(w_idx) = Inf;
        end
    end

end
