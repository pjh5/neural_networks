% DEPRECATED
%
% This file is slated for deletion. Problems with Matlab's approached to object
% oriented programming (in particular, always passing by value) leads to
% annoying difficulty with this approach, and to hacky, non-intuitive code.
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function train_BDH(bdh, x_in, params, step)

    % Developer Notes
    % All matrices are always row based. Data is stored in rows. Prototypes
    %   are stored in rows.
    % data_ prefix implies distances / differences measured in data space
    % diffs_ prefix is for differences (vectors), and dists_ prefix is for
    %   distances (scalars)

    % Set learning rate based on the timestep
    attenuating_learn_rate = max([0.1, params.learning_initial_value * ...
                                    exp(-step / params.learning_decay_time)]);

    % Choose the winning prototype (index of)
    % The winner is the prototype "closest" to the input vector.
    [w_idx, data_dists, data_diffs] = som_winner(bdh.prototypes, x_in);
    winning_prototype = bdh.prototypes(w_idx, :);

    %% Compute Neighborhood Function

    % nbr_sigma is the "width" of the Gaussian neighborhood (which
    % should decrease with timestep).
    % nbr_func is the value of the neighborhood function for every
    % prototype (nPrototypes x 1)
    nbr_sigma = max([params.neighborhood_minimum_width, ...
                    2 * (params.neighborhood_initial_width * ...
                            exp(-step / params.neighborhood_minimum_width))^2]);
    manhattan_dists = sum(abs(...
                bsxfun(@minus, params.grid_indexes, params.of_grid(w_idx))...
                ), 2);
    nbr_func = exp(manhattan_dists .^ 2 / -nbr_sigma^2);

    %% Compute learning rate
    % For the BDH, this is a function of inverse frequency and inverse
    % Voronoi polyhedra size, which is estimated from the distance to the
    % winning prototype.
    bdh_learn_rate = bdh.frequencies(w_idx) * ...
                        data_dists(w_idx) ^ params.estimated_dimensionality;
    learn_rate = attenuating_learn_rate * ...
                                min(0.9, bdh_learn_rate ^ -bdh.magnification);

    % Update prototypes
    bdh.prototypes = bdh.prototypes + learn_rate * ...
                                         bsxfun(@times, nbr_func, data_diffs);

    % Update the frequency of the winning prototypes
    bdh.frequencies(w_idx) = 0;
    bdh.frequencies = bdh.frequencies + ones(params.N_som, 1);
end
