% plot_som_vectors
%
% Plots the weight vectors of the prototypes arranged in the lattice space.
% This is useful when the data space has dimension greater than 3. This can
% help visualize the prototypes and allows comparison between adjacent (in
% lattice space) prototypes.
%
% ARGUMENTS:
%   SOM -- An nP x K matrix of nP prototypes, each of dimension K, arranged
%          into nP rows. These should exist in the data space (should not be
%          scaled into unit space).
%
%   dim -- A vector [width, height] giving the dimensions of the SOM lattice.
%
% RETURNS:
%   A figure with the plot on it.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function fig_handle = plot_som_vectors(SOM, dim)
    fig_handle = figure;
    hold on;
    
    %% Parameters
    margin = 0.05; % distance between side of box and start of vectors
    
    % These are calculated, don't change these calculations
    [N, K] = size(SOM);
    offset = 0.5 - margin; % distance between side of box and vector
    length = offset*2;     % side length of box that encloses vector
    sd2 = dim(2);          % height of prototypes
    
    %% Index mapping functions
    % Maps a prototype index into its (x,y) lattice coordinates
    %   the ugly math is due to matlab's one-indexing
    of_grid = @(i) [ceil(i / sd2), mod(i, sd2) + sd2*(mod(i, sd2)==0)];
    
    %% Draw every box and a vector in every box
    
    % Normalize the prototypes so that they fit in the little boxes
    [f_scale, ~] = scaling_functions(SOM, offset, true);
    SOM = f_scale(SOM);
    
    % Used as the 'x' vector for each prototype
    %   divides the prototypes box (length) into K parts
    %   (1:K) - 1 to correct for one indexing
    %   K - 1 for having an extra "tentpost", e.g. (0, 0.5, 1) for /2
    xv = ((1:K) - 1) / (K-1) * length - offset;
    
    % Loop through every prototype and plot them
    for i=1:N
        loc = of_grid(i);
        x = loc(1);
        y = loc(2);
        
        % Plot the box
        % Around half of the lines are repeated, but it doesn't matter
        % TODO since this is so slow, this probably does matter
        plot([x-.5 x-.5], [y+.5 y-.5], 'k') % left
        plot([x-.5 x+.5], [y+.5 y+.5], 'k') % top
        plot([x+.5 x+.5], [y+.5 y-.5], 'k') % right
        plot([x-.5 x+.5], [y-.5 y-.5], 'k') % bottom
        
        % Plot the vector
        proto = SOM(i,:);
        plot(xv + x, proto + y, 'b');
    end
    
    %% Axes and titles
    xlabel('Prototype Column', 'fontsize', 16);
    ylabel('Prototype Row', 'fontsize', 16);
end
