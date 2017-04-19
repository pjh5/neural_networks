% forward_prop
%
% Propogates the matrix X through the 2-layer MLP represented by the matrices
% W_1 and W_2.
%
% PARAMETERS:
% x             The M x K input matrix of M input vectors each of dimension K
%
% W_1           The N x (K+1) weight matrix mapping from inputs to the first
%               hidden layer
%
% W_2           The Ky x (N+1) weight matrix mapping from the hidden layer to
%               the output layer
%
% dropout_p     (OPTIONAL) The hyperparameter for dropout, this is the
%               percentage of neurons that will be cleared on every layer.
%
% Written by Jesse Hellemn at Rice University, 2016
function [fNET2, fNET1, d2, d1] = forward_prop(X, W_1, W_2, dropout_p)

    % Developer notes
    %
    % This method is written to handle a matrix of data vectors instaed of just
    % 1 vector, even though this is slightly less efficient, so that this same
    % method can be reused for bulk recall for an entire set of vectors (a test
    % set, for example). Note that bulk learning was removed from the train_mlp
    % implementation, because online learning is better.
    %
    % Dropout is implemented separately for efficiency reasons, so that no
    % unnecessary calls to rand are made for non-dropout propogations.

    [M,~] = size(X);

    %% No dropout
    if nargin == 3 || 0 == dropout_p

        % Input to hidden layer
        % [M, K+1] * [K+1, N] => [M, N]
        fNET1 = tanh([ones(M,1) , X] * W_1');

        % Hidden layer to output
        % [M, N+1] * [N+1, Ky] => [M, Ky]
        fNET2 = tanh([ones(M,1) , fNET1] * W_2');

        % No dropout
        d1 = ones(M, size(W_1, 1));
        d2 = ones(M, size(W_2, 1));

    %% Dropout
    else

        % Input to hidden layer
        % [M, K+1] * [K+1, N] => [M, N]
        d1 = (rand(M, size(W_1, 1)) > dropout_p) / dropout_p;
        fNET1 = tanh(d1 .* ([ones(M,1) , X] * W_1'));

        % Hidden layer to output
        % [M, N+1] * [N+1, Ky] => [M, Ky]
        d2 = (rand(M, size(W_2, 1)) > dropout_p) / dropout_p;
        fNET2 = tanh(d2 .* ([ones(M,1) , fNET1] * W_2'));
    end
end
