% plot_layered_prototypes
%
% Equivalent of generic_plot_function for the layered_bdh. The layered_bdh has
% special visualizations. See overlay_layered_prototypes.m and
% plot_layer_trajectories.m for more detailed information.
% 
%
% ARGUMENTS
% X         An N x K matrix of N data points, each of dimension K. Each row is a
%           data vector. K must be 2 or 3. This should not be scaled.
%
% Y         N x 1 vector of class labels, where the ith element is the class
%           label (>= 1) of the ith data vector (ith row of x).
%
% suite     A cell array where each cell is a struct with a field .prototypes.
%           These are the prototypes that will be overlayed above the data.
%           These prototypes are expected to be scaled, as f_unscale will be
%           called on them before plotting. If they are already unscaled, then
%           pass in f_unscale to do nothing.
%
% dims      A vector [width, height] giving the dimensions of each of the BDH
%           lattices (they must all be the same).
%
% f_unscale     A function to unscale the prototypes back into the data space
%
% title_prefix  A string, the prefix to be append to all of the titles.
%
% overlay       Boolean value, whether or not to make an overlay figure.
%
% trajectories  Boolean value, whether or not to plot the trajectories.
%
%
% RETURNS
% Makes a figure if overlay and another if trajectories, but doesn't actually
% return the function handles.
%
% Implemented by Jesse Hellemn at Rice University, Fall 2016
function plot_layered_prototypes(X, Y, suite, dims, f_unscale, ...
                                        title_prefix, overlay, trajectories)

    % Overlay all the prototypes
    if overlay
        [fig_over, leg_plots, leg_str] = overlay_layered_prototypes(...
                                                X, Y, suite, dims, f_unscale);
        title({title_prefix 'Overlayed BDH Layers'}, 'fontsize', 16)
        xlabel('x', 'fontsize', 16)
        ylabel('y', 'fontsize', 16)
        zlabel('z', 'fontsize', 16) % no-op for 2 dims
        legend(leg_plots, leg_str, 'Location', 'northwest')
    end

    % Plot trajectories of every prototype
    if trajectories
        [fig_traj, leg_plots, leg_str] = plot_layer_trajectories(...
                                                X, Y, suite, dims, f_unscale);
        title({title_prefix, 'Prototype Trajectories'}, 'fontsize', 16)
        xlabel('x', 'fontsize', 16)
        ylabel('y', 'fontsize', 16)
        zlabel('z', 'fontsize', 16) % no-op for 2 dims
        legend(leg_plots, leg_str, 'Location', 'northwest')
    end

end
