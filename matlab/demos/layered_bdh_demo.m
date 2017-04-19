% layered_bdh_demo
%
% Example file of how to setup and run a layered_bdh
%
%
% PARAMETERS
% magnification_factors     The magnification factors to train the BDHs with.
%
% effective_dimension   A guess of the "true" dimensionality of the underlying
%                       data manifold.
%
% title_prefix          A string to be used as the prefix of titles on plots.
%
% f_get_data            A function to return an X and Y (data and labels).
%
%
% RETURNS
% suite                 The output of layered_bdh.m
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [suite] = layered_bdh_demo(magnification_factors, ...
                                effective_dimension, title_prefix, f_get_data)

    % Get the data passed in from the user
    [X, Y] = f_get_data();

    %% Parameters
    params = struct();
    dims = [1, 1] * 10;
    params.lattice_dimensions = dims;

    % Iteration parameters
    MAX_ITERS = 20000;
    params.n_max_iterations = MAX_ITERS;
    params.n_ksom_iterations = 3000;

    % Neighborhood parameters
    params.neighborhood_initial_width = max(dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
        1000 / log(params.neighborhood_initial_width);

    % Learning rate parameters
    params.learning_initial_value = 0.1;
    params.learning_decay_time = round(MAX_ITERS);

    % Only BDH parameter is the magnification factor
    params.magnification_factors = magnification_factors;
    params.estimated_dimensionality = effective_dimension;

    % Plotting stuff
    params.f_should_monitor = @(step) step == MAX_ITERS;
    params.f_monitoring = @(suite, step, f_unscale) ...
            plot_layered_prototypes(X, Y, suite, dims, f_unscale, ...
                                                    title_prefix, true, true);

    % Run the SOM
    [suite] = layered_bdh(X, params);

end
