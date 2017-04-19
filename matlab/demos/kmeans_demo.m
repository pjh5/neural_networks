% kmeans_demo
%
% Creates pretty kmeans plots for FCPS datasets for given values of k. In
% particular, the FCPS datasets were designed to be difficult for clustering,
% and this serves as a good demo for how kmeans is not a good clustering
% algorithm. The datasets 'Wingnut', 'Lsun', 'Chainlink', and 'Atom' are
% particularly recommended as examples when kmeans cannot do a good job. For
% an example of when kmeans does (sometimes) work, try 'Hepta', though on
% some runs kmeans will be suboptimal here as well.
%
% FCPS stands for the Fundamental Clustering Problems Suite, first created for
%    Ultsch, A.: Clustering with SOM: U*C, In Proc. Workshop on Self-Organizing
%    Maps, Paris, France, (2005), pp. 75-82
% The datasets are free for download from the Philipps University website at:
% https://www.uni-marburg.de/fb12/arbeitsgruppe
%
%
% ARGUMENTS:
%   dataset -- The name of the FCPS datasets. E.g. 'Atom', or 'Lsun'
%
%   kVals -- An array of all values of k (number of clusters) to run kmeans
%            on. E.g. [1 3 5] or [4]. A plot of a single kmeans run will be
%            created for every value in kvals. You can put duplicates here
%            too if you want to run a value twice (e.g. [1, 5, 5, 5, 7], but
%            the code runs fast enough that you can just call the entire
%            function again.
%
% SIDE EFFECTS:
%   This will produce size(kVals) + 1 plots. One for every value of k, and
%   another plotting mean intra-cluster variance from k=1 to max(kVals).
%
% Written by Jesse Hellemn at Rice University, Fall 2016
function kmeans_demo(kVals, dataset, X, Y)

    % Read Dataset
    if nargin < 3
        [X, Y] = read_fcps(dataset);
    end
    nClasses = max(Y);

    % Parameters for the visualization
    % Marker sizes are the square of the size you want. Default is 6^2
    markerSize = 6 ^ 2;
    centroidSize = (3 ^ 2) * markerSize;

    % Define the colors to use
    % In the many plots that follow, we use markers for the true labels and
    % colors for the cluster labels. The only colors built in are the bright
    % monochromatic and dichromatic colors, since the original author is
    % colorblind. You can additional colors by adding rows to the array below.
    % colors = {'r', 'b', 'y', 'w', 'm', 'c', 'g'};
    colors = [ 1 0 0; 0 0 1; 1 1 0; 1 1 1; 1 0 1; 0 1 1; 0 1 0];
    markers = {'^', 'o', 'p', 's', 'd', 'h', 'v', '<', '>'};
    
    % Complain if too many clusters to color
    maxK = max(kVals);
    nColors = size(colors, 1);
    if maxK > nColors
        disp(['There are ' num2str(maxK) ' clusters but only ' ...
            num2str(nColors) ' colors. Colors will be reused.'])
    end
    colorOf = @(idxs) colors(mod(idxs + nColors - 1, nColors) + 1, :);

    % Use the right plot function based on dimensions of input dataspace
    if size(X, 2) == 2
        f_scat = @(A, mark, mSize) scatter(A(:,1), A(:,2), mSize, ...
                                colorOf(A(:,3)), 'filled', 'Marker', mark);
    elseif size(X, 2) == 3
        f_scat = @(A, mark, mSize) scatter3(A(:,1), A(:,2), A(:,3), ...
                        mSize, colorOf(A(:,4)), 'filled', 'Marker', mark);
    else
        error('kmeans_demo only supports 2 or 3 dimensional data.');
    end

    %% Plots of k-means for several values of k
    for k = kVals
        [clusterOf, centroids] = kmeans(X, k);
        cX = [X clusterOf];

        % Plot the dataset
        data_fig = figure;
        set(gca, 'Color', 'black');
        title([dataset ' with ' num2str(k) ' Clusters'], 'fontsize', 16)
        xlabel('x', 'interpreter', 'latex', 'fontsize', 16)
        ylabel('y', 'interpreter', 'latex', 'fontsize', 16)
        zlabel('z', 'interpreter', 'latex', 'fontsize', 16) % No-op for 2D
        hold on;

        % Plot each true class with a different marker
        for tc = 1:nClasses
            f_scat(cX(Y == tc, :), markers{tc}, markerSize);
        end

        % Plot the centroids
        f_scat([centroids [1:k]'], markers{nClasses + 1}, centroidSize);
    end

    %% Intra Cluster Variance by K
    
    % Calculate intracluster variance for a range of K
    variances = zeros(1, maxK);
    for k = 1:maxK
        [clusterOf, centroids] = kmeans(X, k);

        % Calculate distances to centroids
        quant_errors = l2_norm(X - centroids(clusterOf, :));

        % Calculate mean squared distances by centroid
        variances(k) = 0;
        for cluster = 1:k
            variances(k) = variances(k) + ...
                                var(quant_errors(clusterOf == cluster));
        end
        variances(k) = variances(k) / k;
    end

    % Plot the MSE for a range of k
    figure;
    plot(1:maxK, variances, '-*b', 'linewidth', 3)
    title('Mean Intra-Cluster Variance', 'fontsize', 16)
    xlabel('Number of Clusters (k)', 'fontsize', 16)
    ylabel('Mean Intra-Cluster Variance', 'fontsize', 16)
end


