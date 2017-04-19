% som_winner
%
% Finds the winner for the given input vector. Currently hardcoded to use
% Euclidean distance. 
%
% ARGUMENTS:
% prototypes    An N x K matrix of N prototypes, each of dimension K, arranged
%               into rows.
%
% x_in          A 1 x K input vector, of which to find the closest prototype to.
%
% biases        (OPTIONAL) an N x 1 vector of biases, which will be subtracted
%               from the distances before the winner is chosen. These are used
%               in CSOM.
%
% RETURNS:
% w_idx         The index (row number) of the winning prototype
% data_dists    The L2 norm of data_diffs
% data_diffs    An N x K matrix of (x_in - prototypes)
%
% Written by Jesse Hellemn at Rice University, 2016
function [w_idx, data_dists, data_diffs] = som_winner(prototypes, x_in, biases)

    % Calculate Euclidean differences to all of the prototypes
    data_diffs = bsxfun(@minus, x_in, prototypes);

    % Distance is just the L2 norm of the differences
    data_dists = l2_norm(data_diffs);

    % Biases are optional
    if nargin < 3
        [~, w_idx] = min(data_dists);
    else
        [~, w_idx] = min(data_dists - biases);
    end
end
