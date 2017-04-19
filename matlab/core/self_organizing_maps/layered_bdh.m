% layered_bdh
%
% Co-trains a suite of BDHs. Each BDH has its own magnification factor, but all
% other parameters are the same across them. Each BDH starts with the same
% initial prototype vectors, and are trained with the same order of input
% vectors.
%
% PARAMETERS
% Xraw          The N x K input matrix of N data vectors each of dimension K,
%               arranged into rows. Scaling is not needed, as scaling is done
%               internally.
%
% params        Apart from one extra field, 'magnification_factors', this has
%               the same parameters as BDH. See BDH.m for more details.
%               params.magnification_factors is a vector of magnification
%               factors. One BDH will be trained per magnification factor.
%
%
% RETURNS
% suite     A cell array, with one cell per trained BDH. suite{i} is a struct
%           with two fields:
%           .prototypes     The trained prototypes of the BDH
%           .magnification  The magnification factor of the BDH  
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function [suite] = layered_bdh(Xraw, params)

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

    % Suite parameters
    magnification_factors = params.magnification_factors;
    X_eff_dim = params.estimated_dimensionality;
    n_layers = numel(magnification_factors);

    % Scale data
    % This allows us to visualize the data uniformally. 
    [f_scale, f_unscale] = scaling_functions(Xraw, 1, true);
    X = f_scale(Xraw);

    % Maps a prototype index into its (x,y) lattice coordinates
    %   the ugly math is due to matlab's one-indexing
    lattice_dim = params.lattice_dimensions;
    N_som = prod(lattice_dim);
    latt_height = lattice_dim(2);
    of_grid = @(i) [ceil(i / latt_height), ...
                mod(i, latt_height) + latt_height*(mod(i, latt_height)==0)];
    grid_indexes = of_grid((1:N_som)');

    %% Create every BDH layer

    % Every BDH starts with the same initial weights
    SOM_protos = 2 * rand(N_som, size(Xraw, 2)) - 1;

    % Store frequencies in a 2D matrix
    frequencies = ones(n_layers, N_som);

    % Make every layer
    suite = {};
    for i = 1:n_layers
        suite{i} = struct();
        suite{i}.prototypes = SOM_protos;
        suite{i}.magnification = magnification_factors(i);
    end

    %% Trains with normal KSOM first
    % Note that since this is in a separate loop, the learning rate and
    % neighborhood function will be reset after this. All this should
    % accomplish is a topological ordering of the map
    for step=1:params.n_ksom_iterations
        x_in = X(randi(size(X, 1)),:);
        learn_rate = params.learning_initial_value * ...
                                    exp(-step / params.learning_decay_time);
        nbr_sigma = 2 * (params.neighborhood_initial_width * ...
                            exp(-step / params.neighborhood_decay_time))^2;

        % Train every bdh separately, but with the same input vctor
        for i_bdh = 1:n_layers
            [w_idx, ~, data_diffs] = ...
                                    som_winner(suite{i_bdh}.prototypes, x_in);
            nbr_func = ...
                exp(sum(...
                    abs(bsxfun(@minus, grid_indexes, of_grid(w_idx))) , 2 ...
                       ) .^ 2 / -nbr_sigma ^ 2);

            %% Update prototypes
            suite{i_bdh}.prototypes = suite{i_bdh}.prototypes + learn_rate * ...
                                         bsxfun(@times, nbr_func, data_diffs);

            % Update the frequency of the winning prototypes
            frequencies(i_bdh, w_idx) = 0;
            frequencies(i_bdh, :) = frequencies(i_bdh, :) + ones(1, N_som);
        end
    end


    %% Trains with the BDH algorithm
    % Note that since step starts at 1 again, the learn rate and neighborhood
    % function are reset. This is intentional.
    for step=1:params.n_max_iterations

        % Set learning rate based on the timestep
        attenuating_learn_rate = max([0.01, params.learning_initial_value * ...
                                    exp(-step / params.learning_decay_time)]);

        % nbr_sigma is the "width" of the Gaussian neighborhood (which should
        % decrease with timestep).
        nbr_sigma = max([params.neighborhood_minimum_width, ...
                        2 * (params.neighborhood_initial_width * ...
                            exp(-step / params.neighborhood_decay_time))^2]);

        % Pick a random input and use the same one for every layer
        x_in = X(randi(size(X, 1)),:);

        % Train every bdh separately, but with the same input vctor
        for i_bdh = 1:n_layers

            % Choose the winning prototype (index of)
            % The winner is the prototype "closest" to the input vector.
            [w_idx, data_dists, data_diffs] = ...
                                    som_winner(suite{i_bdh}.prototypes, x_in);

            %% Compute Neighborhood Function

            % nbr_func is the value of the neighborhood function for every
            % prototype (nPrototypes x 1)
            manhattan_dists = sum(abs(...
                    bsxfun(@minus, grid_indexes, of_grid(w_idx))...
                    ), 2);
            nbr_func = exp(manhattan_dists .^ 2 / -nbr_sigma ^ 2);

            %% Compute learning rate
            % For the BDH, this is a function of inverse frequency and inverse
            % Voronoi polyhedra size, which is estimated from the distance to
            % the winning prototype.
            bdh_learn_rate = frequencies(i_bdh, w_idx) * ...
                            data_dists(w_idx) ^ params.estimated_dimensionality;

            % To avoid instability (and to try and avoid flipping of the BDHs
            % away from each other while ordering topographically), keep normal
            % KSOM learning rates until the lattice is topographically ordered
            learn_rate = attenuating_learn_rate * ...
                        min(0.9, bdh_learn_rate ^ -suite{i_bdh}.magnification);

            %% Update prototypes
            suite{i_bdh}.prototypes = suite{i_bdh}.prototypes + learn_rate * ...
                                         bsxfun(@times, nbr_func, data_diffs);

            % Update the frequency of the winning prototypes
            frequencies(i_bdh, w_idx) = 0;
            frequencies(i_bdh, :) = frequencies(i_bdh, :) + ones(1, N_som);
        end

        %% Monitor progress
        if params.f_should_monitor(step)
            params.f_monitoring(suite, step, f_unscale);
        end
    end

    % Unscale the prototypes before they're returned
    for i_bdh = 1:n_layers
        suite{i_bdh}.prototypes = f_unscale(suite{i_bdh}.prototypes);
    end
end
