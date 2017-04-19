% plot_2d_som
%
% Plots a SOM with a 2D lattice in either 2D or 3D data space.
%
% ARGUMENTS:
%   fig -- A figure object created by <code>figure</code> on which to plot
%          the prototypes
%   SOM -- The prototypes of the SOM to plot, organized in an M x K matrix
%          where the ith row is the Kth dimensional weight vector of the ith
%          prototype. 
%   dims - A 2 element vector [x, y] where x is the width and y the height of
%          the SOM lattice.
%   varagin -- Keyword arguments to be passed along to the plot functions.
%
% SIDE EFFECTS:
%   Plots directly onto fig
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function plot_obj = plot_lattice_in_data_space(fig, SOM, dims, varargin)
    figure(fig);
    hold on

    % This code hacks through drawing the 2d lattice to represent the SOM
    % by drawing two "plots". One zigzags through the lattice horizontally,
    % and the other vertically. It's a little hard to write cleanly, but it
    % works. 

    % We start with a simple index grid
    % 1 4 7
    % 2 5 8
    % 3 6 9
    grid_idxs = reshape(1:prod(dims), dims);

    %% Vertical Zigzag

    % Makes a bit mask like
    % 1 0 1
    % 1 0 1
    % 1 0 1
    mask = repmat([ones(dims(1),1) , zeros(dims(1),1)], 1, ceil(dims(2)/2));
    mask = mask(:,1:dims(2));
    
    % Makes alternating vertical indices like
    % 1 6 7
    % 2 5 8
    % 3 4 9
    vert = (grid_idxs.*mask + grid_idxs(dims(1):-1:1, :).*(~mask));
    vert = vert(:);

    %% Horizontal Zigzag

    % Makes a bit mask like
    % 1 1 1
    % 0 0 0
    % 1 1 1
    mask = repmat([ones(1,dims(2)) ; zeros(1,dims(2))], ceil(dims(1)/2), 1);
    mask = mask(1:dims(1),:);

    % Makes alternating horizontal indices something like
    % 1 8 3
    % 6 5 4
    % 9 2 7
    % I'm not actually too sure about this though
    horiz = (grid_idxs.*mask + grid_idxs(:, dims(2):-1:1).*(~mask))';
    horiz = horiz(:);

    %% Call to plot

    % Plot 2 zigzags for fishnet connection
    % Only the first plot object is returned, but this should be fine for
    % legend purposes, since the style and color will be the same for both
    % plots.
    if size(SOM, 2) == 2
        plot_obj = plot(SOM(horiz,1), SOM(horiz,2), varargin{:});
        plot(SOM(vert,1), SOM(vert,2), varargin{:})
    elseif size(SOM, 2) == 3
        plot_obj = plot3(SOM(horiz,1), SOM(horiz,2), SOM(horiz,3), varargin{:});
        plot3(SOM(vert,1), SOM(vert,2), SOM(vert,3), varargin{:})
    end

end
