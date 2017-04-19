% overlay_layered_prototypes
%
% Plots the data with true classes colored, and then every lattice in suite,
% with colors evenly spaced between black and white.
%
%
% ARGUMENTS
% X         An N x K matrix of N data points, each of dimension K. Each row is a
%           data vector. K must be 2 or 3. This should not be scaled.
%
% Y         N x 1 vector of class labels, where the ith element is the class
%           label (>= 1) of the ith data vector (ith row of x).
%
% suite     A cell array where each cell is a struct with a field .prototypes.
%           These are the prototypes that will be overlayed above the data.
%           These prototypes are expected to be scaled, as f_unscale will be
%           called on them before plotting. If they are already unscaled, then
%           pass in f_unscale to do nothing.
%
% dims      A vector [width, height] giving the dimensions of each of the BDH
%           lattices (they must all be the same).
%
%
% f_unscale     A function to unscale the prototypes back into the data space
% 
% 
% RETURNS
% A figure with the requested plot on it.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function [fig, legend_plots, legend_str] = overlay_layered_prototypes(...
                                                X, Y, suite, dims, f_unscale)

    %% Plotting Parameters
    linewidth = 3;
    prototypeSize = 5;
    shadow_size = 3;

    % Useful functions
    n_layers = numel(suite);
    f_color = @(which_layer) [1 1 1] * (which_layer - 1) / n_layers;

    % Overlay all the prototypes

    % Plot the data
    fig = figure;
    plot_labelled_data(fig, X, Y)
    hold on

    % Build the legend string in the right order
    legend_str = {};
    legend_idx = 1;

    % Overlay all the prototype layers
    for i_bdh = n_layers:-1:1
        legend_plots(legend_idx) = plot_lattice_in_data_space(fig, ...
                            f_unscale(suite{i_bdh}.prototypes), dims, ...
                            '-s', ...
                            'Linewidth', linewidth + (i_bdh==1), ...
                            'MarkerSize', prototypeSize, ...
                            'MarkerFaceColor', f_color(i_bdh), ...
                            'Color', f_color(i_bdh));

        legend_str{legend_idx} = num2str(suite{i_bdh}.magnification);
        legend_idx = legend_idx + 1;
    end
end
