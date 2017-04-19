function fig_handle = contour_plot(SOM_dim, wins, step)

    %% Mesh and count of won inputs
    [x, y] = meshgrid(1:SOM_dim(2), 1:SOM_dim(1));
    
    %% Contour
    fig_handle = figure;
    contour(x,y,reshape(sum(wins,2), SOM_dim), 5, 'ShowText', 'on');
    
    %% Axes and labels
    set(gca, 'XTick', 1:10);
    set(gca, 'YTick', 1:10);
    title(['Contour plot of Prototypes Win Counts at Step ' num2str(step)], 'fontsize', 16);
    xlabel('Prototype Column', 'fontsize', 16);
    ylabel('Prototype Row', 'fontsize', 16);
    pause(1);
end