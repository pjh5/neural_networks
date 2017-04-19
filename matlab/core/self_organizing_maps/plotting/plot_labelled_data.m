% plot_labelled_data
%
% Plots the data, all with classes colored. This only works on 2 or 3
% dimensional data space.
%
% ARGUMENTS:
%   X -- An N x K matrix of N data points, each of dimension K. Each row is a
%        data vector. K must be 2 or 3. This should not be scaled.
%
%   Y -- N x 1 vector of class labels, where the ith element is the class
%        label (>= 1) of the ith data vector (ith row of x).
%
% RETURNS:
%   A figure with the requested plot on it.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function plot_labelled_data(fig, X, Y)

    % Parameters for the visualization
    % The marker sizes passed into the matlab plot functions are interpreted as
    % the area of the marker. It's easier to think in terms of radius and
    % square the result.
    markerSize = 4 ^ 2;

    % Markers to preferably use. These markers can be filled, which are much
    % easier to see in 3D plots than non-filled markers (such as *, +, or .).
    % This plot function only uses circle and star, but the entier cell array
    % is here for convenient reference.
    markers = {'^', 'o', 'p', 's', 'd', 'h', 'v', '<', '>'};

    % Map the classes to colors This sample uses a different color for every
    % class. The only colors built in are the bright monochromatic and
    % dichromatic colors, since the original author is colorblind. You can
    % additional colors by adding rows to the array below. 
    % Black is reserved for prototypes
    %  colors = {'r', 'b', 'm', 'y', 'c', 'g'}
    colors = [ 1 0 0; 0 0 1; 1 0 1; 1 1 0; 0 1 1; 0 1 0];

    % If Y is in 1-in-C encoding, change it to integer class labels
    if size(Y, 2) > 1
        [~, Y] = max(Y, [], 2);
    end

    % Complain if there are too many classes to color
    neededColors = max(Y) + 1; % + 1 for the lattice prototypes
    nColors = size(colors, 1);
    if neededColors > nColors
        disp(['There are ' num2str(neededColors) ' clusters but only ' ...
            num2str(nColors) ' colors. Colors will be reused.'])
    end

    % Create a function handle to map labels to colors. This function will map
    % label 1 to the 1st color, label 2 to the 2nd color, etc. Colors will wrap
    % around if there are more needed colors than specified colors.
    colorOf = @(idxs) colors(mod(idxs + nColors - 1, nColors) + 1, :);

    % Use the right plot function based on dimensions of input dataspace
    % X is the data, Y the labels, mark the marker to use (like square or
    % star), and mSize the size of the marker. 'filled' parameter makes the
    % marker filled in, and so won't work on non-fillable markers like * or +
    if size(X, 2) == 2
        f_scat = @(X, Y, mark, mSize) scatter(X(:,1), X(:,2), mSize, ...
                                colorOf(Y), 'filled', 'Marker', mark);
    elseif size(X, 2) == 3
        f_scat = @(X, Y, mark, mSize) scatter3(X(:,1), X(:,2), X(:,3), ...
                        mSize, colorOf(Y), 'filled', 'Marker', mark);
    else
        error('This plot function only handles 2 and 3 dimensional data.')
    end

    % Set up the figure
    figure(fig)
    xlabel('x', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('y', 'interpreter', 'latex', 'fontsize', 16)
    zlabel('z', 'interpreter', 'latex', 'fontsize', 16) % No-op for 2D

    % Plot the dataset
    f_scat(X, Y, markers{2}, markerSize)

end
