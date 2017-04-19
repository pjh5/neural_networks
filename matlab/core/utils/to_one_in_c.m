% to_one_in_c
%
% Converts N x 1 vector of integer class labels to a N x C matrix M, where C
% is the number of classes and M_{i,j} is 1 if and only if data point i has
% class label j. This will no-op if given a matrix of 1-in-C encoded class
% labels, so it should be safe to call this multiple times on the same Y.
%
% ARGUMENTS:
%   vector_labels -- An N x 1 vector of integers, where the ith element is the
%                    class label of the ith data point. If vector_labels is
%                    actually N x K for any K > 1, then this function will
%                    assume that vector_labels is already in a 1-in-C encoding
%                    and will no-op, returning the original vector_labels.
%
%   num_labels    -- The number of classes. The output matrix will have
%                    dimension N x num_labels. This is passed in as a
%                    parameter both for efficiency reasons, and in case
%                    the input labels do not happen to have the highest
%                    labels.
%
% RETURNS:
%   M -- An N x num_labels matrix, where every row has exactly one 1 and every
%        other element is a 0. Every row corresponds to a datapoint.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function M = to_one_in_c(vector_labels, num_labels)

    % No-op if vector_labels is already 1 in C
    if size(vector_labels, 2) > 1
        M = vector_labels;
        return
    end

    % First makes a matrix that looks like
    % 1 2 ... C
    % .
    % .
    % .
    % 1 2 ... C
    % And then the == operator turns it into the wanted array
    M = repmat(1:num_labels, size(vector_labels, 1), 1) == vector_labels;
end
