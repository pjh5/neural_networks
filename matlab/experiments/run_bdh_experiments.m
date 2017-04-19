
% Runs some experiments on BDH alpha and its effects


%% Experiments on the learning rate adjustment
%{

    % This is the function that BDH uses to adjust the learning rate
    f_adjust_lr = @(mag, k, lr_n0, N_som) ...
                        lr_n0 * (sqrt(2)/2 * N_som^(1/k - 1/2)) .^ (k * mag);

    % Plot of the learning rate adjustment for different values of k (X
    % dimensinoality) against magnification factor. This uses an initial
    % specified learning rate of 0.1, and 100 prototypes.
    f_adjust_mk = @(mag, k) f_adjust_lr(mag, k, 0.1, 100);

    % Kvalues to plot for
    kvals  = [  2,   3,   4,   5,   7,  10,  20];
    colors = {'r', 'y', 'g', 'b', 'm', 'c', 'k'};
    legend_str = {['k=' num2str(kvals(1))]};

    % Support of the plot
    xv = -2:0.1:2;

    % Plot one semilogy first to set axis variables before the hold on
    figure
    semilogy(xv, f_adjust_mk(xv, kvals(1)), colors{1})
    hold on

    % Plot the rest of the values of k
    for i = 2:numel(kvals)
        semilogy(xv, f_adjust_mk(xv, kvals(i)), colors{i}, 'linewidth', 2)
        legend_str{i} = ['k=' num2str(kvals(i))];
    end

    % Label the plot
    title({'Adjusted Learning Rate', 'Against Magnification Factor'}, ...
                                                                'fontsize', 16)
    xlabel('Magnification Factor', 'fontsize', 16)
    ylabel('Adjusted Initial Learning Rate', 'fontsize', 16)
    legend(legend_str)
%}

%% Experiments on distribution of BDH learning rates, with and without the bdh
% adjustments. This section plots the 5-number summary (min, IQR, max) of the
% BDH learning rates across several different magnification factors. This is
% dependent on the data, since it is dependent on k, and on the learning
% process of the BDH itself.
    mag_factors = -3 : 0.3 : 3;

    % Parameters to use for all datasets, with or without initial learning rate
    % adjustments
    params = struct();
    MAX_ITERS = 10000;
    lattice_dims = [1, 1] * 10;
    params.lattice_dimensions = lattice_dims;
    params.n_max_iterations = MAX_ITERS;
    params.neighborhood_initial_width = max(lattice_dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
                                1000 / log(params.neighborhood_initial_width);
    params.learning_initial_value = 0.03;
    params.learning_decay_time = round(MAX_ITERS / 2);
    params.f_should_monitor = @(i) i == MAX_ITERS;

    % Atom dataset
    name = 'Atom';
    [X, Y] = read_fcps(name);
    params.estimated_dimensionality = size(X, 2);

    % Monitoring turned off, since past experiments show that the learning is
    % alright. To turn it back on, use the commented out function below this
    make_f_monitoring = @(adjust, mag) @(SOM, step) 1;
    %make_f_monitoring = @(adjust, mag) @(SOM, step) generic_plot_function(...
    %    X, Y, SOM, lattice_dims, ...
    %    [name ' ' repmat('No ', adjust) 'Adjustment, Alpha=' num2str(mag)], ...
    %    step, step == MAX_ITERS, 0, 0);

    [fig_noadj, fig_adj, fig_cap] = measure_learning_rates([0], ...
                            name, X, mag_factors, params, make_f_monitoring);

%% Histogram of learning rates for a typical run. This is just to show the huge
% variety in order of magnitude across learning rates.
%{
    n_bars = 25;

    params = struct();
    MAX_ITERS = 10000;
    lattice_dims = [1, 1] * 10;
    params.lattice_dimensions = lattice_dims;
    params.n_max_iterations = MAX_ITERS;
    params.neighborhood_initial_width = max(lattice_dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
                                1000 / log(params.neighborhood_initial_width);
    params.learning_initial_value = 0.03;
    params.adjust_initial_learning_rate = false;
    params.learning_decay_time = round(MAX_ITERS / 2);
    params.f_should_monitor = @(i) i == MAX_ITERS;
    params.f_monitoring = @(SOM, step) 1; % No monitoring

    % Atom dataset
    name = 'Atom';
    [X, Y] = read_fcps(name);
    params.estimated_dimensionality = size(X, 2);

    % Run a BDH, ignoring everything but learning rates
    params.magnification_factor = 0.5;
    [~, debugs] = BDH(X, params);

    % Plot a histogram of the uncapped learning rates
    % This is complicated in order to get a nice log x-axis.
    fig_unadj_hist = figure;
    unadj_hist = histogram(log(debugs.learning_rates), n_bars);
    xtickformat('10^%g')
    title('Distribution of Uncapped Learning Rates', 'fontsize', 16)
    xlabel('Uncapped Learning Rate (Log Scale)', 'fontsize', 16)
    ylabel('Count', 'fontsize', 16)

    % Plot a histogram of capped learning rates
    fig_adj_hist = figure;
    adj_hist = histogram(min(0.9, debugs.learning_rates), n_bars);
    set(gca, 'yscale', 'log')
    title('Distribution of Capped Learning Rates', 'fontsize', 16)
    xlabel('Capped Learning Rate', 'fontsize', 16)
    ylabel('Count (Log Scale)', 'fontsize', 16)
%}

