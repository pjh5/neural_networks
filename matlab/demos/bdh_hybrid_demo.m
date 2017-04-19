% bdh_hybrid_demo
%
% Example file of how to setup and run a BDH-MLP hybrid. Right now, these
% hybrids will only work on classification problems (the BDH part of the hybrid
% forces vector quantization that probably isn't good for non-classification
% tasks anyways).
%
% PARAMETERS
%
% estimated_dimensionality      A guess of the "true" dimension of the
%                               underlying data manifold. This can be a float.
%                               Whether this really affects the result of
%                               learning the BDH or not isn't well known. But
%                               it probably doesn't hurt, especially if you
%                               expect this to be much less than the dimension
%                               of the input vectors.
%
% X, Y              The input N x K data matrix of N data vectors each of
%                   dimension K, along with the corresponding labels (as
%                   integers)
%
% X_test, Y_test    (Optional) Similar matrix and vector as X and Y, but to
%                   test on. If these are not passed in, then X and Y will be
%                   divided into training and test sets by an 80/20 split.
%
%
% RETURNS
%
% suite             The output of train_bdh_hybrid.m. See there for details
%
% Y_all_hat         A N x NUM_MAGNIFICATION_FACTORS matrix, where
%                   Y_all_hat[i,j] is the jth BDH-MLP hybrid's predicted class
%                   of input vector i.
%
% Written by Jesse Hellemn at Rice University, 2016
function [suite, Y_all_hat] = bdh_hybrid_demo(...
                                estimated_dimensionality, X, Y, X_test, Y_test)

    %% Parameters
    magnification_factors = [1.0, 0.7, 0, -0.7, -1.0];
    %magnification_factors = .7 * ones(1, numel(magnification_factors));
    k_weights = [.6 .2 .2];

    %% Adapt the passed in datasets to training/testing data
        % If only passed one dataset, split it into a train and test sets
        if nargin <= 4
            PERCENT_TRAINING = 0.7;

            % Divide the dataset into a training and test set
            total_size = size(X, 1);
            n_train = round(PERCENT_TRAINING * total_size);
            rand_idxs = randperm(total_size);
            X = X(rand_idxs, :);
            Y = Y(rand_idxs, :);

            X_train = X(1:n_train, :);
            Y_train = Y(1:n_train, :);
            X_test  = X(n_train + 1:end, :);
            Y_test  = Y(n_train + 1:end, :);

        % Otherwise just rename the datasets
        else
            X_train = X;
            Y_train = Y;
        end

        % Classification requires 1-in-C encoding
        n_classes = max(max(Y));
        Y_train = to_one_in_c(Y_train, n_classes);
        Y_test = to_one_in_c(Y_test, n_classes);

        % Finally package data sets into a struct to pass to the MLP
        data = struct();
        data.train        = X_train;
        data.train_labels = Y_train;
        data.test         = X_test;
        data.test_labels  = Y_test;

    %% MLP Parameters
        mlp_params = struct();
        mlp_params.learning_rate        = 0.1;
        mlp_params.momentum             = 0.9;
        mlp_params.dropout_percentage   = 0;
        mlp_params.tolerance            = 0.01;
        mlp_params.N_hidden             = 7;
        mlp_params.max_iterations       = 10000;
        mlp_params.radius_initial_weights = 0.1;

        % Log file
        log_file = fopen(['logs/hybridmlp_' ...
            num2str(mlp_params.N_hidden) '_' ...
            num2str(mlp_params.learning_rate) '_' ...
            num2str(mlp_params.momentum) '_' ...
            num2str(mlp_params.max_iterations) '.txt'], 'w');
        mlp_params.f_log = @(str) fprintf(log_file, str);

        % Monitoring parameters
        monitoring_freq = 100;
        mlp_params.f_should_monitor = @(step) mod(step, monitoring_freq) == 0;
        mlp_params.n_monitoring_samples = ...
                        round(mlp_params.max_iterations / monitoring_freq) + 1;

        % Error function
        mlp_params.f_stopping_error = @(output, desired) ...
                                1 - classification_accuracy(output, desired);

    %% BDH Parameters
        bdh_params = struct();
        bdh_params.lattice_dimensions = [1, 1] * 10;
        bdh_params.n_max_iterations = 30000;
        bhd_params.n_ksom_iterations = 4000;
        bdh_params.estimated_dimensionality = estimated_dimensionality;

        % Neighborhood parameters
        bdh_params.neighborhood_minimum_width = 1;
        bdh_params.neighborhood_initial_width = ...
                                        max(bdh_params.lattice_dimensions) / 2;
        bdh_params.neighborhood_decay_time = ...
                            1000 / log(bdh_params.neighborhood_initial_width);

        % Learning rate parameters
        bdh_params.learning_initial_value = 0.1;
        bdh_params.learning_decay_time = round(bdh_params.n_max_iterations /2);

        % No monitoring. Monitoring is built into the hybrid
        bdh_params.adjust_initial_learning_rate = false;
        bdh_params.f_should_monitor = @(i) false;
        bdh_params.f_monitoring = @(SOM, step) 1;

    %% Train the hybrid
    suite = train_bdh_hybrid(data, magnification_factors, k_weights, ...
                                                    bdh_params, mlp_params);

    %% Collect all propogated outputs for experiments
    n_bdh = numel(magnification_factors);

    % Combine all the data
    X_all = [X_train ; X_test];
    Y_all = [Y_train ; Y_test];

    % Propogate through the BDHs
    X_all_bdh = zeros(...
                size(X_all, 1), prod(bdh_params.lattice_dimensions), n_bdh);
    for i=1:n_bdh
        X_all_bdh(:, :, i) = top_k_bmus(suite{i}.bdh, X_all, k_weights);
    end

    % Propogate through the MLPs
    Y_all_hat = zeros(size(X_all, 1), n_bdh);
    for i=1:n_bdh
        [~, Y_all_hat(:,i)] = max(...
            forward_prop(X_all_bdh(:, :, i), suite{i}.W_1, suite{i}.W_2), ...
            [], 2);
    end

    % Variance of classifier labels
    prediction_variance = var(Y_all_hat, 0, 2);
    fprintf('%f%% of data vectors had non-unanimous classifications\n', ...
                                sum(prediction_variance > 0) / size(X_all, 1))

    % Rounded mean accuracy
    prediction_mean = to_one_in_c(round(mean(Y_all_hat, 2)), n_classes);
    accuracy_mean = classification_accuracy(Y_all, prediction_mean);
    fprintf('When averaging the outputs of all classifiers:\n')
    fprintf('\tMisclassifiaction rate:  %f%%\n', 1 - accuracy_mean)
    disp(confusion_matrix(Y_all, prediction_mean))

    % Mode accuracy
    prediction_mode = to_one_in_c(mode(Y_all_hat, 2), n_classes);
    accuracy_mode = classification_accuracy(Y_all, prediction_mode);
    fprintf('When taking the majority vote of all classifiers:\n')
    fprintf('\tMisclassifiaction rate:  %f%%\n', 1 - accuracy_mode)
    disp(confusion_matrix(Y_all, prediction_mode))
    

end
