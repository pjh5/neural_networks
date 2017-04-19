% BDH
%
% Trains a BDH SOM on Xraw, and returns the prototypes. The
% parameters for the ViSOM are passed in the parameters struct, along with
% callback functions for monitoring.
%
% The algorithm implemented here was introduced by Bauer, Der, and Hermann in
% "Controlling the Magnification Factor of Self-Organizing Feature Maps", 1996.
% See that paper for motivation, reasoning, and further details of this
% algorithm.
%
% ARGUMENTS:
% Xraw      An N x K numeric matrix with N data points, each with K
%           dimensions. Note that every row is a datapoint.
%
% params    Struct that contains the following fields:
%   .magnification
%           The desired magnification factor (theoretically). For 2 dimensional
%           uncorrelated data, Ritter and Schulten showed that this should be
%           1.5 * [desired_magnification] - 1. For correlated data or higher
%           dimensions, there is no good theory to say what this should be.
%
%   .estimated_dimensionality
%           An estimate (possibly a float) of the "true" dimensionality of the
%           underlying data manifold. The learning rate used during learning is
%           calculated using this value, which is used to estimate the mean
%           volume of the receptive fields, which is used as an approximation
%           of the prototype weight vectors' local density in data space. A
%           good estimate may be important to learning, especially for very
%           high dimensional data, but there isn't much theory to say how
%           important this is or what its effect is.
%
%   .adjust_initial_learning_rate
%           A boolean. If set, then this implementation will adjust the passed
%           in learning rate to try an account for the effect of lattice
%           dimensions on the average magnitude of the learning rate. In
%           effect, experimentation on 2 and 3 dimensional data shows that this
%           doesn't do much. This is not in the original BDH paper.
%
%   .lattice_dimensions
%           Numeric vector specifying the dimensions of the SOM lattice.
%           For example, [10, 10] will specify a SOM with 100 prototypes
%           arranged in a square rectangular lattice.
%
%   .n_max_iterations
%           The maximum number of training iterations to carry out
%
%   .neighborhood_initial_width
%           The initial width of the Gaussian neighborhood function. This
%           is approximately equal to the radius around the winning
%           prototype (in lattice space) of neighboring prototypes that
%           will be significantly updated.
%
%   .neighborhood_minimum_width
%           The minimum width that the neighborhood function should use. In
%           most cases, this shosuld be 1, corresponding to always updating
%           immediately neighboring prototypes.
%
%   .neighborhood_decay_time
%           About how long it takes for the neighborhood function to decay
%           to the minimum width
%
% f_should_monitor  A user specified function that should take the
%                   timestep and return either true or false. If true, then
%                   f_monitoring will be called. This function will be called
%                   on every timestep.  f_should_monitor(timestep)
%
% f_monitoring      A user specified monitoring function that will be
%                   called on EVERY iteration. f_monitoring is passed the
%                   prototypes (scaled back to original data space), the
%                   current timestep.  f_monitoring(prototypes, timestep)
%
% RETURN VALUE:
% SOM_protos    A nP x K numeric matrix containing nP prototypes, each
%               of dimension K. Each row stores a prototype.  Technically,
%               these are just the weight vectors in the data space that are
%               associated with each prototype. The "coordinates" of the
%               prototypes in the lattice space are stored implicitly by their
%               row-index in the matrix.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [SOM_protos, debugs] = BDH(Xraw, params)

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
    lattice_dim = params.lattice_dimensions;
    N_som = prod(lattice_dim);
    max_iters = params.n_max_iterations;
    nbr_sigma_0 = params.neighborhood_initial_width;
    nbr_min_width = params.neighborhood_minimum_width;
    nbr_t1 =  params.neighborhood_decay_time;

    % Monitoring functions
    f_should_monitor = params.f_should_monitor;
    f_monitoring = params.f_monitoring;

    % Controls the magnification factor
    magnification = params.magnification_factor;
    X_eff_dim = params.estimated_dimensionality;

    % Learning rate for the BDH
    % For BDH, the learning rate is actually determined at every iteration
    % based on estimates of data density around the winning prototype. We keep
    % the learning rate passed in on the struct as the same learning rate used
    % in other SOM variants (the same scale, usually around 0.1-0.01). But the
    % dynamically decided learning rate from the BDH varies with the frequency
    % and Voronoi cell volume of the winning prototype, which in turn varies
    % with the size of the SOM and the width of the data.
    lr_n0 = params.learning_initial_value;
    lr_t2 = params.learning_decay_time;
   
    % Adjust the learning rate to account for the potential difference in
    % magnitude of the BDH learning rates. This is an optional check so that we
    % can run experiments on this adjustment.
    if params.adjust_initial_learning_rate
        learn_rate_0 = lr_n0 * (sqrt(2)/2 * N_som^(1/X_eff_dim - 1/2)) ...
                                                ^ (X_eff_dim * magnification);
        disp(sprintf('The learning rate was changed to %d', learn_rate_0))
    else
        learn_rate_0 = lr_n0;
    end

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
    SOM_protos = 2*rand(N_som, size(Xraw, 2)) - 1;
    SOM_idxs = (1:N_som)';
    SOM_grid_idxs = of_grid(SOM_idxs);
    frequencies = ones(N_som, 1);

    % Debugging stuff
    d_freqs = zeros(max_iters, 1);
    d_vols = zeros(max_iters, 1);
    d_learns = zeros(max_iters, 1);

    %% Trains the SOM
    for i=1:max_iters

        % Set learning rate based on the timestep
        attenuating_learn_rate = max([learn_rate_0 * exp(-i / lr_t2), 0.1]);

        % Pick a random input and find its closest prototype
        % Randomly pick an input vector at every timestep. 
        % This samples with replacement, which should be the desired
        % behavior. With replacement is much more efficient, and it doesn't
        % really matter if repeats are encountered.
        x_in = X(randi(size(Xraw, 1)),:);

        % Choose the winning prototype (index of)
        % The winner is the prototype "closest" to the input vector.
        [w_idx, data_dists, data_diffs] = som_winner(SOM_protos, x_in);
        winning_prototype = SOM_protos(w_idx, :);

        %% Compute Neighborhood Function
        % nbr_sigma is the "width" of the Gaussian neighborhood (which
        % should decrease with timestep).
        % nbr_func is the value of the neighborhood function for every
        % prototype (nPrototypes x 1)
        nbr_sigma = max([2*(nbr_sigma_0 * exp(-i / nbr_t1))^2, nbr_min_width]);
        lattice_diffs = bsxfun(@minus, SOM_grid_idxs, of_grid(w_idx));
        manhattan_dists = sum(abs(lattice_diffs), 2);
        nbr_func = exp(manhattan_dists .^ 2 / -nbr_sigma^2);

        %% Compute learning rate
        % For the BDH, this is a function of inverse frequency and inverse
        % Voronoi polyhedra size, which is estimated from the distance to the
        % winning prototype.
        bdh_learn_rate = frequencies(w_idx) * data_dists(w_idx)^X_eff_dim;
        d_learns(i) = bdh_learn_rate ^ -magnification;
        final_learn_rate = attenuating_learn_rate * ...
                                    min(0.9, bdh_learn_rate ^ -magnification);

        % Debugging
        d_freqs(i) = frequencies(w_idx) ^ -1;
        d_vols(i) = data_dists(w_idx) ^ -X_eff_dim;

        % Update prototypes
        SOM_protos = SOM_protos + final_learn_rate * ...
                                bsxfun(@times, nbr_func, data_diffs);

        % Update the frequency of the winning prototypes
        frequencies(w_idx) = 0;
        frequencies = frequencies + ones(N_som, 1);

        %% Monitor progress
        if f_should_monitor(i)
            f_monitoring(f_unscale(SOM_protos), i);
        end
    end

    % Unscale the prototypes before they're returned
    SOM_protos = f_unscale(SOM_protos);

    % Debugging
    debugs = struct();
    debugs.frequencies = d_freqs;
    debugs.volumes = d_vols;
    debugs.learning_rates = d_learns;
end
