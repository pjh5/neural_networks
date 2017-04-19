% CSOM
%
% Trains a CSOM (Conscience SOM) on Xraw, and returns the prototypes. The
% parameters for the CSOM are passed in the parameters struct, along with
% callback functions for monitoring.
%
% This can only handle rectangular 2D lattices.
%
% The algorithms implemented in this file are due to <>
%
%
% ARGUMENTS:
%     Xraw - An N x K numeric matrix with N data points, each with K
%            dimensions. Note that every row is a datapoint.
%
%     params - Struct that contains the following fields:
%     .lattice_dimensions
%           Numeric vector specifying the dimensions of the SOM lattice.
%           For example, [10, 10] will specify a SOM with 100 prototypes
%           arranged in a square rectangular lattice.
%
%     .n_max_iterations
%           The maximum number of training iterations to carry out
%
%     .learning_initial_value
%           The initial value of the learning rate. Experimentally, 0.1
%           seems to work pretty well.
%
%     .learning_decay_time
%           About how many iterations it takes for the learning rate to
%           decay to small values. Past this timestep, the SOM will change
%           very very slowly.
%
%     .csom_neighborhood_width
%           The exact width of the neighborhood, in prototypes. This defaults
%           to 1, as was used in the original paper.
%
%     .csom_bias_weight
%           The weight that the winning frequency of a prototype will affect
%           its winning, relative to its distance from the input prototype.
%           When this number is 0, the algorithm reduces to KSOM with a
%           constant fixed width neighborhood function.
%
%     .csom_bias_decay
%           This is the rate at which "winning frequency" of a prototype
%           decays. Specifically, this amount of the winning frequency will
%           decay every iteration, except for the prototype that wins (which
%           will have 1 added to its winning frequency). This should really be
%           a function of how many data vectors and how many prototypes there
%           are, as the number of times a prototype wins is propotional to
%           (number of data points / number of prototypes). That said, the CSOM
%           algorithm appears pretty robust to variations in this value, though
%           it should be kept relatively small.
%
%     f_should_monitor - A user specified function that should take the
%                        timestep and return either true or false. If true,
%                        then f_monitoring will be called. This function will
%                        be called on every timestep.
%                           f_should_monitor(timestep)
%
%     f_monitoring - A user specified monitoring function that will be
%                    called on EVERY iteration. f_monitoring is passed the
%                    prototypes (scaled back to original data space), the
%                    current timestep.
%                        f_monitoring(prototypes, timestep)
%
% RETURN VALUE:
%     SOM_protos - A nP x K numeric matrix containing nP prototypes, each
%                  of dimension K. Each row stores a prototype.
%                  Technically, these are just the weight vectors in the
%                  data space that are associated with each prototype. The
%                  "coordinates" of the prototypes in the lattice space are
%                  stored implicitly by their row-index in the matrix.
%
% Implementation by Jesse Hellemn at Rice University, Fall 2016
function SOM_protos = CSOM(Xraw, params)

    % Developer Notes
    % All matrices are always row based. Data is stored in rows. Prototypes
    %   are stored in rows.
    % data_ prefix implies distances / differences measured in data space
    % diffs_ prefix is for differences (vectors), and dists_ prefix is for
    %   distances (scalars)
    % All internal computations are done on scaled data (scaled between -1
    %   and 1. But the prototype matrix passed to the monitoring function
    %   and prototypes returned from this function are all unscaled back to
    %   the original data space.

    %% Data Input and Parameters
    [N_in, X_dim] = size(Xraw);
    lattice_dim = params.lattice_dimensions;
    N_som = prod(lattice_dim);
    max_iters = params.n_max_iterations;
    lr_n0 = params.learning_initial_value;
    lr_t2 = params.learning_decay_time;

    % CSOM neighborhood width, by default is 1
    % In the original paper, it was 1
    if isfield(params, 'csom_neighborhood_width')
        neighborhood_width = params.csom_neighborhood_width;
    else
        neighborhood_width = 1;
    end

    % Relative weight of winning frequency relative to distance from input
    bias_constant = params.csom_bias_weight;

    % This should be a function of data size, since with more data points wins
    % will be spaced more out and decay will happen more between wins
    bias_decay = params.csom_bias_decay;

    % Monitoring functions
    f_should_monitor = params.f_should_monitor;
    f_monitoring = params.f_monitoring;

    % Scale data
    % This allows us to visualize the data uniformally. 
    [f_scale, f_unscale] = scaling_functions(Xraw, 1, true);
    X = f_scale(Xraw);

    % Maps a prototype index into its (x,y) lattice coordinates
    %   the ugly math is due to matlab's one-indexing
    latt_height = lattice_dim(2);
    of_grid = @(i) [ceil(i / latt_height), ...
        mod(i, latt_height) + latt_height*(mod(i, latt_height)==0)];

    %% Create the SOM prototypes
    SOM_protos = 2*rand(N_som, X_dim) - 1;
    SOM_idxs = (1:N_som)';
    SOM_grid_idxs = of_grid(SOM_idxs);

    % Allocate win frequencies
    % For the CSOM, prototypes keep track of their recent win frequency. These
    % are not the biases, but biases are calculated using these values
    proto_freqs = zeros(N_som, 1);

    %% Trains the SOM
    for step=1:max_iters

        % Set learning rate based on the timestep
        learn_rate = max([lr_n0 * exp(-step / lr_t2), 0.1]);

        % Pick a random input and find its closest prototype
        % Randomly pick an input vector at every timestep. 
        % This samples with replacement, which should be the desired
        % behavior. With replacement is much more efficient, and it doesn't
        % really matter if repeats are encountered.
        x_in = X(randi(N_in),:);

        % Calculate the biases
        biases = bias_constant * (1/N_som - proto_freqs);

        % Choose the winning prototype (index of)
        % The winner is the prototype "closest" to the input vector.
        [w_idx, ~, data_diffs] = som_winner(SOM_protos, x_in, biases);

        % Update the winning frequencies of the prototypes
        proto_freqs = (1 - bias_decay) * proto_freqs;
        proto_freqs(w_idx) = proto_freqs(w_idx) + bias_decay;

        % Compute Neighborhood Function
        % For the CSOM, this is just the adjacent vectors to the winner
        lattice_diffs = bsxfun(@minus, SOM_grid_idxs, of_grid(w_idx));
        nbr_func = sum(abs(lattice_diffs), 2) <= neighborhood_width;

        % Update prototypes
        SOM_protos = SOM_protos + learn_rate * ...
                                bsxfun(@times, nbr_func, data_diffs);

        %% Plot
        if f_should_monitor(step)
            f_monitoring(f_unscale(SOM_protos), step);
        end
    end

    % Unscale the prototypes before they're returned
    SOM_protos = f_unscale(SOM_protos);
end
