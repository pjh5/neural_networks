% som_demo
%
% This is a sample of how to use the ViSOM interface, and some sample (good)
% parameters for it.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function SOM = som_demo(dataset, X, Y)

    % Read in data, either from FCPS or passed from the user
    if nargin < 1
        error('This function needs either the name of a FCPS dataset or a matrix of data X and its labels Y')
    elseif nargin < 2
        [Xd, Yd] = read_fcps(dataset);
    else
        Xd = X;
        Yd = Y;
    end

    %% Parameters
    % All the parameters on the struct must be included, and cannot be renamed.
    params = struct();
    lattice_dims = [15, 15];
    params.lattice_dimensions = lattice_dims;

    % Iteration parameters
    % free_iterations is actually ignored right now TODO
    % n_ksom_iterations is the number of "normal" KSOM before specialized ViSOM or
    % CSOM kick in.
    MAX_ITERS = 20000;
    params.n_free_iterations = 2000;
    params.n_max_iterations = MAX_ITERS;
    params.ksom_iterations = 2 * params.n_max_iterations; % Normal KSOM

    % Neighborhood parameters
    params.neighborhood_initial_width = max(lattice_dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
        1000 / log(params.neighborhood_initial_width);

    % Learning rate parameters
    params.learning_initial_value = 0.1;
    params.learning_decay_time = params.n_max_iterations;

    % ViSOM parameters must be present if calling ViSOM
    params.visom_resolution = 0.5;

    % CSOM parameters must be present if calling CSOM
    params.csom_bias_weight = 10;
    params.csom_bias_decay = 0.00001;

    %% Training / Monitoring parameters

    % Define when to plot
    % This will plot the prototypes in data space at step 1000 and the end,
    % never plot histograms, and plot the modified uMatrix at the end.
    when_protos = @(step) step == 1000 || step == MAX_ITERS;
    when_hists = @(step) false;
    when_uMat = @(step) step == MAX_ITERS;
    params.f_should_monitor = @(i) when_protos(i) || when_hists(i) || when_uMat(i);

    % Define what to plot
    params.f_monitoring = @(SOM, step) generic_plot_function(...
                        Xd, Yd, SOM, lattice_dims, dataset, step, ...
                        when_protos(step), when_hists(step), when_uMat(step));

    % Run the SOM
    % This will create the plots if f_plot creates plots.
    SOM = ViSOM(Xd, params);
end
