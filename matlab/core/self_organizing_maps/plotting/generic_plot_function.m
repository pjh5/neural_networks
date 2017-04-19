% generic_plot_function
%
% Wrapper that makes calls to specified plot functions. The last three
% parameters are booleans which tell this which plotting functions to call.
% This can be used as a monitoring function for the ViSOM function.
%
% ARGUMENTS:
%   X -- An N x K matrix of N data vectors, each of dimension K, arranged into
%        N rows.
%
%   Y -- An N x 1 vectors of class labels, where the ith element is the class
%        label of the ith data vector (the ith row of X).
%
%   SOM -- An nP x K matrix of nP prototypes, each of dimension K, arranged
%          into nP rows. 
%
%   dims -- A vector [width, height] giving the dimensions of the SOM lattice.
%
%   title_prefix -- The name of the dataset, used only in the titles of the plots.
%
%   show_protos -- Boolean, whether or not plot the prototypes. If the data
%                  exists in 2 or 3 dimensional space, then the prototypes will
%                  be plotted in data space, superimposed over the data. If the
%                  data space is at least 4 dimensional, then the prototype
%                  vectors are plotted, arranged in the lattice space.
%
%   show_hists -- Boolean, whether or not to plot the histogram of class
%                 labels won by each prototype.
%
%   show_umat -- Boolean, whether or not to plot a modified u-Matrix
%
% RETURNS:
%   Nothing. The figures aren't returned because the number of figures depends
%   on the last three parameters.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function generic_plot_function(X, Y, SOM, dims, title_prefix, step, ...
                                    show_protos, show_hists, show_umat)

    % Compute the number of classes
    % This assumes that the classes are labelled consecutively starting from 1
    nClasses = max(max(Y));

    % Plot the prototypes somehow
    if show_protos

        % 2 or 3 dimensionsl, plot in original data space
        if size(X, 2) < 4
            data_fig = figure;
            plot_labelled_data(data_fig, X, Y)
            hold on
            plot_lattice_in_data_space(data_fig, SOM, dims, '-sk', ...
                    'Linewidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 5)
            title([title_prefix ' ' num2str(step) ' Steps'], 'fontsize', 16)
        else
            data_fig = plot_som_vectors(SOM, dims);
            title([title_prefix ' SOM Vectors ' num2str(step) ' Steps'], ...
                                                        'fontsize', 16)
        end
    end

    % Plot histograms for every prototype
    if show_hists
        hist_fig = figure;
        title(['Prototype Wins for ' title_prefix num2str(step) ' Steps'], ...
                                                            'fontsize', 16)

        [~, win_sums] = full_recall(SOM, X, to_one_in_c(Y, nClasses));
        som_hists(hist_fig, win_sums, dims);
    end

   
    % Plot a modified u-Matrix
    if show_umat
       modu_fig = mod_u_matrix(X, Y, SOM, dims);
       title([title_prefix ' Modified U-Matrix ' num2str(step) ' Steps'], ...
                                                            'fontsize', 16)
    end

end
