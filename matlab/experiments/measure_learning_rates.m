function [fig_no_adj, fig_adj, fig_cap] = measure_learning_rates(adjusts, ...
                            name, X, mag_factors, params, make_f_monitoring)

    fig_no_adj = figure;
    fig_adj = figure;
    percent_capped = zeros(1, numel(mag_factors));
    percent_small  = zeros(1, numel(mag_factors));
    percent_vsmall = zeros(1, numel(mag_factors));

    for adjust = adjusts
        params.adjust_initial_learning_rate = adjust;
        if adjust
            summary_fig = fig_adj;
        else
            summary_fig = fig_no_adj;
        end

        % Keep track of the 5 number summary of BDH learning rates across
        % magnification factors
        lr_0    = zeros(1, numel(mag_factors));
        lr_25   = zeros(1, numel(mag_factors));
        lr_50   = zeros(1, numel(mag_factors));
        lr_75   = zeros(1, numel(mag_factors));
        lr_100  = zeros(1, numel(mag_factors));
        lr_mean = zeros(1, numel(mag_factors));

        % Run the BDH algorithm for each magnification factor
        for i = 1:numel(mag_factors)

            % Run the BDH
            params.magnification_factor = mag_factors(i);
            params.f_monitoring = make_f_monitoring(adjust, mag_factors(i));
            [~, debugs] = BDH(X, params);

            % Calculate the 5 number summary of the learning rates
            sorted_lrs = sort(debugs.learning_rates);
            n_lrs = numel(sorted_lrs);

            lr_0(i)    = sorted_lrs(1);
            lr_25(i)   = sorted_lrs(round(n_lrs / 4));
            lr_50(i)   = sorted_lrs(round(n_lrs / 2));
            lr_75(i)   = sorted_lrs(round(3 * n_lrs / 4));
            lr_100(i)  = sorted_lrs(n_lrs);
            lr_mean(i) = mean(mean(sorted_lrs));

            if ~adjust
                percent_capped(i) = sum(sorted_lrs >= 0.9) / (1.0*n_lrs);
                percent_small(i)  = sum(sorted_lrs < 0.1) / (1.0*n_lrs);
                percent_vsmall(i) = sum(sorted_lrs < 0.01) / (1.0*n_lrs);
            end
        end

        %% Plot the BDH learning rates 5 number summary
        figure(summary_fig)
        semilogy(mag_factors, lr_mean, '-bp')
        hold on
        semilogy(mag_factors, lr_100, '--bv')
        semilogy(mag_factors, lr_75, '-bv')
        semilogy(mag_factors, lr_50, '-bo')
        semilogy(mag_factors, lr_25, '-b^')
        semilogy(mag_factors, lr_0, '--b^')

        % Label the 5 number summary
        if adjust
            title({'Summary of Learning Rates' name}, 'fontsize', 16)
        else
            title({'Summary of Learning Rates' [name ', No Adjustment']}, ...
                                                                'fontsize', 16)
        end
        xlabel('Magnification Factor', 'fontsize', 16)
        ylabel('Learning Rate (Log Scale)', 'fontsize', 16)
        legend({'Mean', 'Max', '75%', 'Median', '25%', 'Min'}, ...
                                                        'Location', 'north')
        grid on
        gcacopy = gca;
        extraxes = axes('Position', get(gcacopy, 'Position'), ...
                    'Color','none', 'XTick',[], 'YAxisLocation','right');
        linkaxes([gcacopy extraxes], 'xy');
    end

    %% Plot the percentage of capped learning rates
    fig_cap = figure;
    plot(mag_factors, percent_capped, 'b', 'linewidth', 2)
    hold on
    plot(mag_factors, percent_small, 'r', 'linewidth', 2)
    plot(mag_factors, percent_vsmall, 'k', 'linewidth', 2)
    title({'Percentage of Learning Rates' 'Pushed to Extremes'}, 'fontsize', 16)
    xlabel('Magnification Factor', 'fontsize', 16)
    ylabel('Percentage of Learning Rates', 'fontsize', 16)
    legend({'Capped at 0.9', '< 0.1', '< 0.01'}, 'Location', 'northwest')

end
