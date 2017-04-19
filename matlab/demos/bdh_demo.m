% bdh_demo
%
% Example file of how to setup and run a BDH.
%
%
% PARAMETERS
% magnification         The magnification factor to train the BDH with.
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
% prototypes    The prototypes of the trained BDH.
%
% debugs        A struct with various helpful debugging stuff
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [prototypes, debugs] = bdh_demo(...
                magnification, effective_dimension, title_prefix, f_get_data)

    % Get the data passed in from the user
    [X, Y] = f_get_data();

    %% Parameters
    params = struct();
    lattice_dims = [1, 1] * 7;
    params.lattice_dimensions = lattice_dims;

    % Iteration parameters
    MAX_ITERS = 100000;
    params.n_max_iterations = MAX_ITERS;

    % Neighborhood parameters
    params.neighborhood_initial_width = max(lattice_dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
        1000 / log(params.neighborhood_initial_width);

    % Learning rate parameters
    params.learning_initial_value = 0.03;
    params.learning_decay_time = round(MAX_ITERS / 2);
    params.adjust_initial_learning_rate = 0;

    % Only BDH parameter is the magnification factor
    params.magnification_factor = magnification;
    params.estimated_dimensionality = effective_dimension;

    %% Training / Monitoring parameters

    % Define when to plot
    % This will plot the prototypes in data space at step 1000 and the end,
    % never plot histograms, and plot the modified uMatrix at the end.
    when_protos = @(step) step == MAX_ITERS;
    when_hists = @(step) false;
    when_uMat = @(step) step == MAX_ITERS;
    params.f_should_monitor = ...
                        @(i) when_protos(i) || when_hists(i) || when_uMat(i);

    % Define what to plot
    params.f_monitoring = @(prototypes, step) generic_plot_function(...
                        X, Y, prototypes, lattice_dims, title_prefix, step, ...
                        when_protos(step), when_hists(step), when_uMat(step));

    % Run the SOM
    % This will create the plots if f_plot creates plots.
    [prototypes, debugs] = BDH(X, params);
    
    % Debug stuff
    lrs = sort(debugs.learning_rates);
    disp(sprintf('Min    : %d', lrs(1)));
    disp(sprintf('25%%    : %d', lrs(round(MAX_ITERS / 4))));
    disp(sprintf('Median : %d', lrs(round(MAX_ITERS / 2))));
    disp(sprintf('Mean   : %d', mean(debugs.learning_rates)));
    disp(sprintf('75%%    : %d', lrs(3 * round(MAX_ITERS / 4))));
    disp(sprintf('Max    : %d', lrs(MAX_ITERS)));

end
