% csom_demo
%
% This is a sample of how to use the CSOM interface, and some sample (good)
% parameters for it.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
close all
all_datasets = {};

%% Thick shell atom
dataset = 'Atom';
noise_var = 2;
proportions = [4, 1];
[Xd, Yd] = read_fcps(dataset, proportions, noise_var);

all_datasets{1} = struct();
all_datasets{1}.dataset = 'Dense Shell Atom';
all_datasets{1}.Xd = Xd;
all_datasets{1}.Yd = Yd;
all_datasets{1}.kvals = [5];

%% Thick nucleus atom
dataset = 'Atom';
noise_var = 2;
proportions = [0, 3];
[Xd, Yd] = read_fcps(dataset, proportions, noise_var);

all_datasets{2} = struct();
all_datasets{2}.dataset = 'Dense Nucleus Atom';
all_datasets{2}.Xd = Xd;
all_datasets{2}.Yd = Yd;
all_datasets{2}.kvals = [5];

%% Image cube
dataset = '8 Class Imagecube';
[Xd, Yd] = read_imagecube();
all_datasets{3} = struct();
all_datasets{3}.dataset = dataset;
all_datasets{3}.Xd = Xd;
all_datasets{3}.Yd = Yd;
all_datasets{3}.kvals = [8];

%% WDBC
formatspec = ...
        '%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
wdbc = readtable(...
        '~/pjh5/neural_nets/data/breast_cancer_wisconsin/wdbc.csv', ...
        'Delimiter', ',', 'Format', formatspec);
Xwdbc = table2array(wdbc(:, 3:end));
Ywdbc = 1 + (table2array(wdbc(:, 2)) == 'B');
all_datasets{4} = struct();
all_datasets{4}.dataset = 'WDBC';
all_datasets{4}.Xd = Xwdbc;
all_datasets{4}.Yd = Ywdbc;
all_datasets{4}.kvals = [5];


%% Reused parameters
params = struct();
lattice_dims = [10, 10];
params.lattice_dimensions = lattice_dims;

% Iteration parameters
% free_iterations is actually ignored right now TODO
% n_ksom_iterations is the number of "normal" KSOM before ViSOM kicks in
MAX_ITERS = 30000;
params.n_free_iterations = 2000;
params.n_max_iterations = MAX_ITERS;
params.ksom_iterations = 2 * params.n_max_iterations; % Normal KSOM

% Neighborhood parameters

% Learning rate parameters
params.learning_initial_value = 0.1;
params.learning_decay_time = params.n_max_iterations;

% ViSOM parameters must be present if calling ViSOM
params.visom_resolution = 0.5;

% Define when to plot
% This will plot the prototypes in data space at step 1000 and the end,
% never plot histograms, and plot the modified uMatrix at the end.
when_protos = @(step) step == MAX_ITERS;
when_hists = @(step) false;
when_uMat = @(step) step == MAX_ITERS;
params.f_should_monitor = @(i) when_protos(i) || when_hists(i) || when_uMat(i);
                
%% Run all experiments on all datasets
for dset = 1:size(all_datasets, 2)-1
    dataset = all_datasets{dset}.dataset;
    Xd = all_datasets{dset}.Xd;
    Yd = all_datasets{dset}.Yd;
    kvals = all_datasets{dset}.kvals;
    
    %% CSOM
    params.neighborhood_initial_width = 1;
    params.neighborhood_minimum_width = 0.5;
    params.neighborhood_decay_time = ...
                            1000 / log(params.neighborhood_initial_width);
    params.csom_bias_weight = 10;
    params.csom_bias_decay = 0.00001;
    params.f_monitoring = @(SOM, step) generic_plot_function(...
                    Xd, Yd, SOM, lattice_dims, [dataset ' Conscience'], step, ...
                    when_protos(step), when_hists(step), when_uMat(step));
    CSOM(Xd, params);
%{
    %% KSOM
    params.neighborhood_initial_width = max(lattice_dims) / 2;
    params.neighborhood_minimum_width = 1;
    params.neighborhood_decay_time = ...
                            1000 / log(params.neighborhood_initial_width);
    params.f_monitoring_function = @(SOM, step) generic_plot_function(...
                    Xd, Yd, SOM, lattice_dims, [dataset ' Kohonen'], step, ...
                    when_protos(step), when_hists(step), when_uMat(step));
    ViSOM(Xd, params, f_should_plot, f_plot);

    %% Kmeans
    if size(Xd, 2) <= 3
        kmeans_demo(kvals, dataset, Xd, Yd);
    end
%}
end
