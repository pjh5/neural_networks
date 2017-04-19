% ViSOM
%
% Trains a ViSOM (Visual SOM) on Xraw, and returns the prototypes. The
% parameters for the ViSOM are passed in the parameters struct, along with
% callback functions for monitoring. This can be parametrized to perform
% normal KSOM, by setting the number of KSOM iterations to at least the
% max iterations, with a very minimal penalty in efficiency.
%
% This can only handle rectangular 2D lattices.
%
% The algorithms implemented in this file are due to Hujun Yin from:
%   Yin, Hujun (2002). Data visualisation and manifold mapping using the
%   ViSOM. Neural Networks, 15, 1005-1016
% This implementation follows the algorithms given in the paper exactly,
% except for the implementation of the resolution parameter. In fact, the
% handling of the resolution parameter here may be incorrect, and wasn't
% actually tested on different dataspaces.
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
%     .n_free_iterations
%           The number of training iterations to carry out before
%           decreasing the neighborhood and learning rate
%
%     .n_max_iterations
%           The maximum number of training iterations to carry out
%
%     .neighborhood_initial_width
%           The initial width of the Gaussian neighborhood function. This
%           is approximately equal to the radius around the winning
%           prototype (in lattice space) of neighboring prototypes that
%           will be significantly updated.
%
%     .neighborhood_minimum_width
%           The minimum width that the neighborhood function should use. In
%           most cases, this shosuld be 1, corresponding to always updating
%           immediately neighboring prototypes.
%
%     .neighborhood_decay_time
%           About how long it takes for the neighborhood function to decay
%           to the minimum width
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
%     .non_visom_iterations
%           The number of timesteps of normal Kohonen SOM learning to use
%           before applying ViSOM. Before this timestep, learning is
%           completely identical to a normal SOM. To disable ViSOM and use
%           normal SOM, simply set this to the maximum number of iterations
%           specified earlier.
%
%     .resolution
%           The resolution of the ViSOM, this seems to correspond to about
%           the distance in data space between adjacent prototypes (ViSOM
%           attempts to enforce about the same distance between adjacent
%           prototypes). With lower resolution, more prototypes might be
%           needed to cover the entire spread of data. With too high of a
%           resolution, the prototypes will be pushed out of the manifold
%           of the data, resulting in many topology violations.
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
% Written by Jesse Hellemn at Rice University, Fall 2016
function SOM_protos = ViSOM(Xraw, params)

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
    free_iters = params.n_free_iterations;
    max_iters = params.n_max_iterations;
    nbr_sigma_0 = params.neighborhood_initial_width;
    nbr_min_width = params.neighborhood_minimum_width;
    nbr_t1 =  params.neighborhood_decay_time;
    lr_n0 = params.learning_initial_value;
    lr_t2 = params.learning_decay_time;
    ksom_iters = params.ksom_iterations;
    resolution = params.visom_resolution;

    % Monitoring functions
    f_monitoring = params.f_monitoring;
    f_should_monitor = params.f_should_monitor;

    % Scale data
    % This allows us to visualize the data uniformally. 
    [f_scale, f_unscale] = scaling_functions(Xraw, 1, true);
    X = f_scale(Xraw);

    % This maps the resolution to the right scale
    % The resolution passed in by the user is sorta proportional and close to
    % the actual distances between prototypes in data space that the ViSOM
    % will attempt to conform to. The resolution used internally must be
    % scaled to handle differences of scale between dataspace and lattice
    % space.
    resolution = resolution * 4 * sqrt(max(var(X))) / min(lattice_dim);
    
    % Maps a prototype index into its (x,y) lattice coordinates
    %   the ugly math is due to matlab's one-indexing
    latt_height = lattice_dim(2);
    of_grid = @(i) [ceil(i / latt_height), ...
        mod(i, latt_height) + latt_height*(mod(i, latt_height)==0)];
    
    %% Create the SOM prototypes
    SOM_protos = 2*rand(N_som, X_dim) - 1;
    SOM_idxs = (1:N_som)';
    SOM_grid_idxs = of_grid(SOM_idxs);

    %% Trains the SOM
    max_iters = max([max_iters, free_iters]);
    for i=1:max_iters
        
        % Set learning rate based on the timestep
        learn_rate = max([lr_n0 * exp(-i / lr_t2), 0.1]);
        
        %% Pick a random input and find its closest prototype
        % Randomly pick an input vector at every timestep. 
        % This samples with replacement, which should be the desired
        % behavior. With replacement is much more efficient, and it doesn't
        % really matter if repeats are encountered.
        x_in = X(randi(N_in),:);
        
        % Choose the winning prototype (index of)
        % The winner is the prototype "closest" to the input vector.
        [w_idx, ~, ~] = som_winner(SOM_protos, x_in);
        winning_prototype = SOM_protos(w_idx, :);
        
        %% Compute Neighborhood Function
        % nbr_sigma is the "width" of the Gaussian neighborhood (which
        % should decrease with timestep).
        % nbr_func is the value of the neighborhood function for every
        % prototype (nPrototypes x 1)
        nbr_sigma = max([2*(nbr_sigma_0 * exp(-i / nbr_t1))^2,   
                         nbr_min_width
                         ]);
        lattice_diffs = bsxfun(@minus, SOM_grid_idxs, of_grid(w_idx));
        manhattan_dists = sum(abs(lattice_diffs), 2);
        nbr_func = exp(manhattan_dists .^ 2 / -nbr_sigma^2);
        
        %% Compute Distances for non-Winning Prototypes
        % For the ViSOM, adjustements are actually made by decomposing the
        % distance between any prototype v and the input vector x (with
        % winning prototype w) into
        %    d(v, x) = d(x, w) + (c * d(v, w))
        % Where d(a,b) is the distance between a and b in the data space,
        % and c is a contraction_modifier, which restricts the contractive
        % fores in the SOM
        data_diff_winner_to_input = x_in - winning_prototype;
        data_diffs_to_winner = bsxfun(@minus, ...
                                        winning_prototype, SOM_protos);

        %% Compute contraction modifiers for the ViSOM
        if i > ksom_iters
            data_dists_to_winner = l2_norm(data_diffs_to_winner);
            contraction_modifiers = (data_dists_to_winner ...
                                        ./ (manhattan_dists * resolution)) - 1;

            % One small quirk - the contraction_modifier for the winner is now
            % infinity, since the lattice_distance from the winner to itself is
            % 0. Since the contraction_modifier isn't really applied to the
            % winner, we just set that particular value to one
            contraction_modifiers(w_idx) = 1;

            contract_diffs = bsxfun(@times, ...
                            contraction_modifiers, data_diffs_to_winner);

        % Not using ViSOM yet, then for normal KSOM there are no
        % contraction modifiers.
        else 
            contract_diffs = data_diffs_to_winner;
        end

        %% Update all prototypes
        % Re-combine the decomposed differences from each prototype to the
        % input vecotr. Remember that d(v,x) = d(w,x) + (c * d(v,w))
        % The learning rate and neighboring function are also factored in
        proto_updates = bsxfun(@plus, ...
                                data_diff_winner_to_input, contract_diffs);
        proto_updates = learn_rate * ...
                                bsxfun(@times, nbr_func, proto_updates);

        % Update prototypes
        SOM_protos = SOM_protos + proto_updates;

        %% Plot
        if f_should_monitor(i)
            f_monitoring(f_unscale(SOM_protos), i);
        end
    end

    % Unscale the prototypes before they're returned
    SOM_protos = f_unscale(SOM_protos);
end
