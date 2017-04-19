% mod_u_matrix
%
% Plots a modified U-matrix for the SOM specified by prototypes and
% lattice_dims. The returned figure will have axis labels, but these can
% always be overwritten.
%
% ARGUMENTS
%   X -- An N x K matrix on N data vectors arranged in rows, each of dimension
%        K. X(i,:) should be the ith data vector.
%
%   Y -- An N x 1 matrix of class labels, where Y(i) is the class label of
%        data vector X(i,:).
%
%   prototypes -- An M x K matrix of M prototypes arranged in rows, each of
%                 which exists in K dimensional space.
%
%   lattice_dims -- A vector [lattice_width, lattice_height]
%
% RETURNS:
%   A figure handle to the figure created.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function fig_handle = mod_u_matrix(X, Y, prototypes, lattice_dims)
    fig_handle = figure;
    hold on;

    % Extract densities from X and Y
    [~, values] = full_recall(prototypes, X, to_one_in_c(Y, max(max(Y))));

    %% Index mapping functions
    sd1 = lattice_dims(1);
    sd2 = lattice_dims(2);

    % Maps a prototype index into its (x,y) lattice coordinates
    %   the ugly math is due to matlab's one-indexing
    of_grid = @(i) [ceil(i / sd2), mod(i, sd2) + sd2*(mod(i, sd2)==0)];

    % Maps a prototype's (x,y) lattice coordinates into its index in the
    % prototypes
    of_arr = @(x, y) (x-1)*sd2 + y;

    %% Fences
    % The fence is the divider between adjacent boxes. The whiter the
    % fence, the larger the distance between the corresponding prototypes.

    % Index groups
    [xb,yb] = meshgrid(1:(sd1), 1:(sd2-1));
    [xu,yu] = meshgrid(1:(sd1), 2:(sd2));
    [xl,yl] = meshgrid(1:(sd1-1), 1:(sd2));
    [xr,yr] = meshgrid(2:(sd1), 1:(sd2));

    % Calculate distances between adjacent prototypes
    f_l2norm = @(V) sqrt(sum(V .^ 2, 2));
    dv = [xl(:) yl(:) f_l2norm(prototypes(of_arr(xr(:), yr(:)),:) - ...
                                    prototypes(of_arr(xl(:), yl(:)),:))];
    dh = [xb(:) yb(:) f_l2norm(prototypes(of_arr(xu(:), yu(:)),:) - ...
                                    prototypes(of_arr(xb(:), yb(:)),:))];

    % Normalize Differences
    maxd = 0.9 * max(max([dv(:,3) ; dh(:,3)]));
    dv(:,3) = dv(:,3) ./ maxd;
    dh(:,3) = dh(:,3) ./ maxd;

    % Draw every vertical line
    for i=1:size(dv, 1)
        x = dv(i,1);
        y = dv(i,2);
        plot([x+.5 x+.5], [y+.5 y-.5], ...
                        'Color', [0 0 0]+min([dv(i,3), 1]), 'linewidth', 5)
    end

    % Draw every horizontal line
    for i=1:size(dh,1)
        x = dh(i,1);
        y = dh(i,2);
        plot([x-.5 x+.5], [y-.5 y-.5], ...
                        'Color', [0 0 0]+min([dh(i,3), 1]), 'linewidth', 5)
    end

    %% Draw every box

    % Adjust counts to be 3 dimensional
    n_counts = size(values, 1);
    counts_dim = size(values, 2);

    % More than 3 dimensions get mapped to red scale (sum of counts)
    if counts_dim == 1 || counts_dim > 3
        values = sum(values, 2);
        values = values / max(values);
        values = [values zeros(n_counts, 1) zeros(n_counts, 1)];

    % 2 dimensions counts get mapped to red and blue
    elseif counts_dim == 2
        values = bsxfun(@rdivide, values, max(values, [], 1));
        values = [values(:, 1) zeros(n_counts, 1) values(:, 2)];

    % 3 dimensions mapped to rgb
    else
        values = bsxfun(@rdivide, values, max(values, [], 1));
    end

    % Plot every box
    for i=1:size(values,1)
        loc = of_grid(i);
        x = loc(1);
        y = loc(2);
        v = [x-.45 y+.45 ; x+.45 y+.45 ; x+.45 y-.45 ; x-.45 y-.45];
        f = [1 2 3 4];
        patch('Faces', f, 'Vertices', v, 'FaceColor', values(i,:));
    end

    %% Axes and titles
    set(gca, 'XTick', 1:sd1);
    set(gca, 'YTick', 1:sd2);
    xlabel('Prototype Column', 'fontsize', 16);
    ylabel('Prototype Row', 'fontsize', 16);
end
