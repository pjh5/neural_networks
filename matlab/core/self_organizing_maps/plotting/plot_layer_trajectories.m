% plot_layer_trajectories
%
% Plots the trajectories of the prototypes across magnification factors. This
% plotting function is specific to layered_bdh.m. The presumption is that suite
% contains several BDHs which have been cotrained with the same starting
% prototypes and the same order of input vectors, and so the prototypes across
% BDHs are linked together by their position in the lattice. This plots the
% first lattice (assumed to be the highest magnification factor), and then
% "trajectories" from each prototype to its location in the next BDH, in
% lighter and lighter shades of grey.
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

% Implemented by Jesse Hellemn at Rice University, Fall 2016
function [fig, legend_plots, legend_str] = plot_layer_trajectories(...
                                                X, Y, suite, dims, f_unscale)

    %% Plotting Parameters
    linewidth = 3;
    prototypeSize = 5;
    shadow_size = 3;

    % Use lighter colors for trajectories 
    n_layers = numel(suite);
    f_color = @(which_layer) [1 1 1] * (which_layer + 1) / (n_layers + 2);

    % Plot the data
    fig = figure;
    plot_labelled_data(fig, X, Y)
    hold on

    % Don't plot the prototypes yet, or the trajectories will be plotted on
    % top of them

    % Collect all of the prototypes into a multidimensional matrix
    [N_protos, X_dim] = size(suite{1}.prototypes);

    % Pick plot function based on dimensionality
    if X_dim == 2
        plot_connection = @(i_bdh, i_proto, prev_p, next_p) ...
                       plot([prev_p(i_proto, 1) next_p(i_proto, 1)], ...
                            [prev_p(i_proto, 2) next_p(i_proto, 2)], ...
                            '-s', ...
                            'Linewidth', linewidth, ...
                            'MarkerSize', shadow_size, ...
                            'Color', f_color(i_bdh) ...
                           );
    elseif X_dim == 3
        plot_connection = @(i_bdh, i_proto, prev_p, next_p) ...
                       plot3([prev_p(i_proto, 1) next_p(i_proto, 1)], ...
                            [prev_p(i_proto, 2) next_p(i_proto, 2)], ...
                            [prev_p(i_proto, 3) next_p(i_proto, 3)], ...
                            '-s', ...
                            'Linewidth', linewidth, ...
                            'MarkerSize', shadow_size, ...
                            'Color', f_color(i_bdh) ...
                           );
    else
        error('This plot function only handles 2 and 3 dimensional data.')
    end

    % Build the legend entries in the right order
    legend_str = {};
    legend_idx = 1;

    % Plot the connection from every prototype to its successor
    % Move in reverse order so that stronger trajectories are on top
    for i_bdh = (n_layers - 1):-1:1
        prev_p = f_unscale(suite{i_bdh}.prototypes);
        next_p = f_unscale(suite{i_bdh + 1}.prototypes);

        for i_proto = 1:N_protos
            legend_plots(legend_idx) = plot_connection(...
                                            i_bdh, i_proto, prev_p, next_p);
        end

        legend_str{legend_idx} = num2str(suite{i_bdh}.magnification);
        legend_idx = legend_idx + 1;
    end

    % Plot the first layer last, to be on top of the trajectories
    legend_plots(legend_idx) = plot_lattice_in_data_space(fig, ...
                            f_unscale(suite{1}.prototypes), dims, ...
                            '-s', ...
                            'Linewidth', linewidth + 1, ...
                            'MarkerSize', prototypeSize, ...
                            'MarkerFaceColor', f_color(-1), ...
                            'Color', f_color(-1));
    legend_str{legend_idx} = num2str(suite{1}.magnification);


end
