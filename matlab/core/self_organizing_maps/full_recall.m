% full_recall
%
% Computes the winner of every data point and collects the won vectors for
% every prototype.
%
% ARGUMENTS:
%   prototypes -- An M x K matrix of M prototypes, each of which is a row
%                 vector that exists in K dimensional space.
%
%   X -- An N x K matrix of data points, each of which is a row vector with
%        dimensionality K.
%
%   Y -- An N x C matrix of class labels, where the ith row is the class
%        label of the ith data point. For this function, it makes most
%        sense for these class labels to be in 1 in C encodings, for
%        example [0 1 0] to represent having class label 2 out of 3
%        possible classes.
%
% RETURNS:
%   winner_idxs -- An N x 1 vector, where the ith element is the index of
%                  the winning prototype for the ith data point.
%
%   won_sums -- An M x C numerical matrix where M is the number of
%               prototypes and C is the dimensionality of the class labels
%               Y. The ith row is obtained by summing the class labels of
%               all the data points that are won by the ith prototype. For
%               this to make sense, the class labels should be vectors,
%               such as the 1 in C encodings, and not just one dimensional
%               integers.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [winner_idxs, won_sums] = full_recall(prototypes, X, Y)
    N_in = size(X,1);

    % Allocate memory to store results
    winner_idxs = zeros(N_in, 1);
    won_sums = zeros(size(prototypes, 1), size(Y, 2));

    % Find the winner for every input data vector
    for i=1:N_in

        % Choose winner
        [w_idx, ~, ~] = som_winner(prototypes, X(i,:));

        % Record winner
        winner_idxs(i) = w_idx;
        won_sums(w_idx,:) = won_sums(w_idx,:) + Y(i,:);
    end
end
