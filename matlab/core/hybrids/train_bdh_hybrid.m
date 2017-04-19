% train_bdh_hybrid
%
% Trains a suite of unconnected BDH-MLP hybrids, with given magnification
% factors. The BDH-MLP hybrids are trained completely separately from each
% other. The prototypes of the BDH and weight vectors of the MLP are returned
% in a cell array, one cell per BDH-MLP hybrid.
%
% This only works for classification problems. In the future it may be adapted
% to handle non-classification problems, but this might not lead to good
% results. The BDH part of the hybrid also leads to vector quantization of the
% original inputs, so real valued outputs (from the hybrid) are unlikely to be
% that accurate.
%
%
% PARAMETERS
%
% data      See the comments for train_mlp.m
%
% magnification_factors     A vector of magnification factors. One BDH-MLP
%                           hybrid will be trained for each magnification
%                           factor in this vector.
%
% k_weights                 A vector of length K of floats between 0 and 1.
%                           This determines the weights of the top K BMUs of
%                           the data after it's mapped through the BDH. See
%                           top_k_bmus.m for more detailed information.
%
% bdh_params                Parameter struct for the BDH. See BDH.m
%
% mlp_params                Parameter struct for the MLP. See train_mlp.m
%
%
% RETURNS
%
% suite     A cell array with the same number of cells as magnification
%           factors. suite{i} is a struct with fields:
%
%           .bdh        The prototypes of the trained BDH. See BDH.m
%
%           .W_1        The weight matrix of the first layer of the MLP, see
%                       train_MLP.m
%
%           .W_2        The weight matrix of the second layer of the MLP, see
%                       train_MLP.m
%
%           .errors     The training and testing errors of training the MLP,
%                       see train_MLP.m
%
% Written by Jesse Hellemn at Rice Universtiy, 2016
function suite = train_bdh_hybrid(data, magnification_factors, k_weights, ...
                                                        bdh_params, mlp_params)

    % Train every BDH hybrid separately in a loop
    suite = {};
    for i = 1:numel(magnification_factors)
        magnification = magnification_factors(i);

        fprintf('Training hybrid for magnification %d\n', magnification)
        bdh_params.magnification_factor = magnification;
        suite{i} = struct();

        % Train a BDH with this magnification factor
        suite{i}.bdh = BDH(data.train, bdh_params);

        % Plot everything to make sure the BDH learned correctly
        % TODO allow this to be specified through parameters somehow
        generic_plot_function(data.train, data.train_labels, suite{i}.bdh, ...
            bdh_params.lattice_dimensions, ...
            ['Alpha=' num2str(magnification) ' '], ...
            bdh_params.n_max_iterations, true, false, false)

        % Map the data through the BDH prototypes
        % NOTE this isn't scaled right now. The MLP requires that the input
        % data is scaled into the range of hyperbolic tangent; This data should
        % already by in the range of hyperbolic tangent (because it will be
        % only the values in k_weights, or 0), but it has not been optimally
        % separated into that range.
        mlp_data = struct();
        mlp_data.train = top_k_bmus(suite{i}.bdh, data.train, k_weights);
        mlp_data.test  = top_k_bmus(suite{i}.bdh, data.test,  k_weights);
        mlp_data.train_labels = data.train_labels;
        mlp_data.test_labels  = data.test_labels;

        % Train the MLP on the output of the BDH
        [W_1, W_2, errors] = train_mlp(mlp_data, mlp_params);
        suite{i}.W_1 = W_1;
        suite{i}.W_2 = W_2;
        suite{i}.errors = errors;

        % Plot the learning curves of the MLP to make sure it learned correctly
        % TODO make this a parameter as well.
        figure
        hold on
        plot(errors.steps, errors.training, 'linewidth', 3)
        plot(errors.steps, errors.testing, 'linewidth', 3)
        hold off

        title({'Learning History' ...
            ['for Magnification: ' num2str(magnification)]}, 'fontsize', 16)
        xlabel('Learning Step', 'fontsize', 16)
        ylabel('Misclassification Rate', 'fontsize', 16)
        legend({'Training Error', 'Test Error'})

        % Display the confusion matrices to give more insight into the learning
        fprintf('\nMagnification = %d, Confusion Matrix:\n', magnification)
        test_confusion_mat = confusion_matrix(...
                mlp_data.train_labels, forward_prop(mlp_data.train, W_1, W_2));
        disp(test_confusion_mat)
        fprintf('\n')
    end
end
