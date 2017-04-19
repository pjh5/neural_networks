% mlp_demo
%
% Trains a MLP on the given data
%
% PARAMETERS:
%   title_prefix -- The name of the dataset, this is only used in plot titles.
%
%   classification -- Boolean (1 or 0) if this is a classification problem or
%                     not. Note that if this is a classification problem, then
%                     the Y should be integers corresponding to class labels
%                     (or 1-in-C encodings of class labels); having floats in Y
%                     and asking for classification will cause problems.
%
%   X -- An N x K matrix of N data vectors, each of dimension K. If the last
%        two parameters X_test and Y_test are present, then this will be used
%        as the training data, otherwise it will be randomly shuffled and split
%        into training and test data.
%
%   Y -- One of two things:
%        - If this is a classification problem, then Y should be the class
%        labels of X, either as a 1-in-C encoding or as integers.
%        - If this is not a classification problem, then this should be a matrix
%        with N rows, and of arbitrary number of columns.
%
%
%
% Written by Jesse Hellemn at Rice University, 2016
function [W_1, W_2, errors] = mlp_demo(title_prefix, classification, ...
                                                        X, Y, X_test, Y_test)

    %% Parametrs
    PERCENT_TRAINING    = 0.7;
    learning_rate       = 0.1;
    momentum            = 0.9;
    dropout_p           = 0;
    tolerance           = 0.001;
    N_hidden            = 20;
    max_iterations      = 40000;

    % Copy parameters to struct
    params = struct();
    params.learning_rate        = learning_rate;
    params.momentum             = momentum;
    params.dropout_percentage   = dropout_p;
    params.tolerance            = tolerance;
    params.N_hidden             = N_hidden;
    params.max_iterations       = max_iterations;

    % Rarely changed parameters
    params.radius_initial_weights = 0.1;
    radius_data = 0.85;

    % Log file
    log_file = fopen(['logs/' title_prefix '_mlp_', num2str(N_hidden), '_', ...
        num2str(learning_rate), '_', num2str(momentum), '_', ...
        num2str(max_iterations), '.txt'], 'w');
    params.f_log = @(str) fprintf(log_file, str);

    % Monitoring parameters
    monitoring_freq = 100;
    params.f_should_monitor = @(step) mod(step, monitoring_freq) == 0;
    params.n_monitoring_samples = round(max_iterations / monitoring_freq) + 1; 

    %% Adapt the passed in datasets
    % If only passed one dataset, split it into a train and test sets
    if nargin < 6

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
    if classification
        n_classes = max(max(Y));
        Y_train = to_one_in_c(Y_train, n_classes);
        Y_test = to_one_in_c(Y_test, n_classes);
    end
        
    % Scale data
    % Data has to be scaled into the range of hyperbolic tangent, which is
    % the non-linearity used in the MLP implementation
    [f_scale_x, ~] = ...
                scaling_functions([X_train ; X_test], radius_data, true);

    X_train = f_scale_x(X_train);
    X_test = f_scale_x(X_test);

    % For not classification problems, the Y must be scaled too
    if ~classification
        [f_scale_y, f_unscale_y] = ...
                scaling_functions([Y_train ; Y_test], radius_data, true);

        Y_train = f_scale_y(Y_train);
        Y_test = f_scale_y(Y_test);
    end

    % Finally package data sets into a struct to pass to the MLP
    data = struct();
    data.train        = X_train;
    data.train_labels = Y_train;
    data.test         = X_test;
    data.test_labels  = Y_test;

    %% Error function
    % This is defined here because it potentially needs the scaling functions
    if classification
        f_err = @(output, desired) 1 - classification_accuracy(output, desired);
    else
        f_err = @(output, desired) ...
            sqrt(mean(mean((f_unscale_y(output) - f_unscale_y(desired)).^2)));
    end
    params.f_stopping_error = f_err;

    %% Train the MLP
    [W_1, W_2, errors] = train_mlp(data, params);
    fclose(log_file);


    %% Display Results

    % Learning curves
    figure
    hold on
    plot(errors.steps, errors.training, 'linewidth', 3)
    plot(errors.steps, errors.testing, 'linewidth', 3)
    hold off

    title([title_prefix 'Learning History'], 'fontsize', 16)
    xlabel('Learning Step', 'fontsize', 16)
    if classification
        ylabel('Classification Accuracy', 'fontsize', 16)
    else
        ylabel('RMSE Error', 'fontsize', 16)
    end
    legend({'Training Error', 'Test Error'})
    print(['images/' title_prefix '_' num2str(N_hidden), '_', ...
        num2str(learning_rate), '_', ...
        num2str(momentum), '.png'], '-dpng');

    % Confusion matrices if classification problem
    if classification
        train_confusion_mat = confusion_matrix(...
                                    Y_train, forward_prop(X_train, W_1, W_2));
        test_confusion_mat = confusion_matrix(...
                                    Y_test, forward_prop(X_test, W_1, W_2));

        % Print confusion matrices
        fprintf('\n\n Training Data Confusion Matrix:\n')
        disp(train_confusion_mat)
        fprintf('\n Testing Data Confusion Matrix:\n')
        disp(test_confusion_mat)

        confusion_file = fopen(['logs/cmat_', num2str(N_hidden), '_', ...
                num2str(learning_rate), '_', num2str(momentum), '.txt'], 'w');
        fprintf(confusion_file, 'Training Data Confusion Matrix:\n')
        print_matrix(confusion_file, train_confusion_mat)
        fprintf(confusion_file, '\nTest Data Confusion Matrix:\n')
        print_matrix(confusion_file, test_confusion_mat)
        fclose(confusion_file);
    end
end
