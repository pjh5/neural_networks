% BPLearn
%
% Trains a 2-layer, feed-forward, multilayer perceptron, trained with
% backpropogation, to learn the mapping between X and Y.
%
% This particular implementation uses hyperbolic tangent as its non-linearity.
% Both momentum and dropout are implemented, but optional (by setting the
% corresponding hyperparameter to 0).
%
% This function always uses a train and a test set (if test sets are not passed
% in, then the passed in datasets are randomly shuffled and split into train
% and test sets); the test sets are only used to monitor learning in learning
% curves.
%
% All inputs to this function (data vectors and their class labels) should be
% scaled to be between -1 and 1, since range of the output of each perceptron
% is between -1 and 1 (because the non-linearity used in hyperbolic tangent).
% Actually, because the slope of the non-linearity attenuates to 0 as it's
% value approaches -1 or 1, it's recommended that all inputs be scaled between
% -0.9 and 0.9, or -0.85 and 0.85 or some similar values. Scaling is not done
% internally, as it is done for all the SOM variants, since the inverse scaling
% functions will probably be needed by client functions. In the case of the
% SOM, the prototypes could be similary scaled back to the data space before
% they are returned, but it is not possible to simply scale the weight matrices
% that define each layer of the MLP.
%
% This does not allow bulk learning, because online learning is better.
%
% PARAMETERS
% data          A struct with the following fields
% .train        The training set, a matrix of vectors with observations
%               arranged in rows
%
% .train_labels The class labels of data.train, in 1-in-C encoding, also
%               arranged in rows
%
% .test and .test_labels, analagous to the .train versions but for the test
%                set. These will not factor into the learning at all, but are
%                only needed for learning curves.
%
% Written by Jesse Hellemn at Rice University, Spring 2016
function [W_1, W_2, errors] = train_mlp(data, params)

    % Unpack Sets
    X_train = data.train;
    Y_train = data.train_labels;
    X_test  = data.test;
    Y_test  = data.test_labels;

    % Initialize variables
    [M, X_dim] = size(X_train);
    [~, Y_dim] = size(Y_train);

    % Unpack learning parameters
    learn_rate = params.learning_rate;
    momentum = params.momentum;
    dropout_p = params.dropout_percentage;
    tolerance = params.tolerance;
    N_hidden = params.N_hidden;
    max_iterations = params.max_iterations;

    % Radius of initial random weights and of data
    radius_initial_weights = params.radius_initial_weights;

    % Estimate on number of monitoring samples that will be taken
    n_monitoring_samples = params.n_monitoring_samples;
    f_should_monitor = params.f_should_monitor;

    % Error functions
    f_stopping_err = params.f_stopping_error;

    % Logging file to use to save data
    f_log = params.f_log;

    %% Initialize matrices and scale data

    % Initialize random weight matrices
    rand_weight = @(m, n) 2 * radius_initial_weights * rand(m, n) - ...
                                                    radius_initial_weights;
    W_1 = rand_weight(N_hidden, X_dim + 1);    % [N  , K+1]
    W_2 = rand_weight(Y_dim, N_hidden + 1);    % [Ky , N+1]
    W_1p = zeros(size(W_1));                   % [N  , K+1]
    W_2p = zeros(size(W_2));                   % [Ky , N+1]

    % Error matrices
    errors = struct();
    errors.steps = zeros(1, n_monitoring_samples);
    errors.training = zeros(1, n_monitoring_samples);
    errors.testing = zeros(1, n_monitoring_samples);
    n_samples = 0;

    %% Loop until error is small enough
    for step = 1:max_iterations

        % Pick a random input
        x_idx = randi(M);
        x = X_train(x_idx, :);
        y = Y_train(x_idx, :);


        %% Forward Propogate

        % [1, K] => [1, Ky]
        % fNET1 = [1, N]
        % fNET2 = [1, Ky]
        [fNET2, fNET1, drop2, drop1] = forward_prop(x, W_1, W_2, dropout_p);


        %% Backward Propogate

        % [1, Ky] - [1, Ky] .* [1, Ky] => [1, Ky]
        delta_2 = drop2 .* (y - fNET2) .* (1 - fNET2.^2);

        % [1, N] .* [1, Ky] * [Ky, N] => [1, N] 
        delta_1 = drop1 .* (1 - fNET1.^2) .* (delta_2 * W_2(:, 2:end));


        %% Update Weights

        % Calculate deltas
        % [N, 1] * [1, K+1] => [N, K+1]
        W_1d = learn_rate * (delta_1' * [1 , x] + momentum * W_1p);

        % [Ky, 1] * [1, N+1] => [Ky, N+1]
        W_2d = learn_rate * (delta_2' * [1 , fNET1] + momentum * W_2p);

        % Weight update
        W_1 = W_1 + W_1d;
        W_2 = W_2 + W_2d;

        % Save deltas for use in momentum of next iteration
        W_1p = W_1d;
        W_2p = W_2d;


        %% Check for convergence and monitor learning periodically
        if step >= max_iterations || f_should_monitor(step)

            % Calculate learning error on train and test sets
            train_err = f_stopping_err(...
                                    forward_prop(X_train, W_1, W_2), Y_train);
            test_err  = f_stopping_err(...
                                    forward_prop(X_test,  W_1, W_2), Y_test);

            % Save errors for later analysis
            n_samples = n_samples + 1;
            errors.training(n_samples) = train_err;
            errors.testing(n_samples) = test_err;
            errors.steps(n_samples) = step;

            % Log
            f_log(sprintf(['Step %10d: training error is' ...
                '%10f, testing error is %10f\n'], step, train_err, test_err));

            % Stop if error is low enough
            if train_err < tolerance
                disp(sprintf('Desired tolerance reached by step %d', step))
                break
            end
        end % Convergence check and monitoring
    end % Training for loop

    % Strip excess zeros off of the error vectors
    errors.training = errors.training(1, 1:n_samples);
    errors.testing  = errors.testing(1, 1:n_samples);
    errors.steps    = errors.steps(1, 1:n_samples);

    % Print very summarized results to stdout
    fprintf('\nTraining stopped at %d steps\n', step);
    fprintf('Final training error: %10f\n', errors.training(n_samples));
    fprintf('Final testing error: %10f\n', errors.testing(n_samples));

    % Log more complete results to the log file
    f_log(sprintf('\n\nTraining stopped at %d steps\n', step));
    f_log(sprintf('Final training error: %10f\n', errors.training(n_samples)));
    f_log(sprintf('Final testing error: %10f\n', errors.testing(n_samples)));
end
