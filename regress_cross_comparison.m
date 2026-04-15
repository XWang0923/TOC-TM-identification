function [t, ret] = regress_cross_comparison(cutoffs, south_ratios, north_ratios)

figure('Visible', 'off');
t = tiledlayout(length(cutoffs), length(cutoffs), 'TileSpacing', 'compact');
ret = nan(length(cutoffs), length(cutoffs), 6);
for k1 = 1:length(cutoffs)
    for k2 = 1:length(cutoffs)
        nexttile;
        plot([0 100], [1 1], 'k--', 'LineWidth', 1, 'Color', [1 1 1] * 0.5); hold on;
        x = [0 5 10 25 50 75 90 95 100];
        y1 = south_ratios(k1, k2, 1:9); y1 = y1(:); y1(y1<=0) = nan;
        scatter(x(y1>=1), y1(y1>=1), [], 'bv', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3); hold on;
        scatter(x(y1<1), y1(y1<1), [], 'rv', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3); hold on;
        y2 = north_ratios(k1, k2, 1:9); y2 = y2(:); y2(y2<=0) = nan;
        scatter(x(y2>=1), y2(y2>=1), [], 'b^', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3); hold on;
        scatter(x(y2<1), y2(y2<1), [], 'r^', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3); hold on;
        left_mean = geomean([y1(1:4); y2(1:4)], [], 'omitmissing');
        right_mean = geomean([y1(6:9); y2(6:9)], [], 'omitmissing');
        ret(k1, k2, :) = [geomean([y1(5:9)'; y2(5:9)'], [], 'omitmissing'), right_mean];
        if left_mean < 1
            plot([0 25], [left_mean left_mean], 'r-', 'LineWidth', 2); hold on;
        else
            plot([0 25], [left_mean left_mean], 'b-', 'LineWidth', 2); hold on;
        end
        if right_mean < 1
            plot([75 100], [right_mean right_mean], 'r-', 'LineWidth', 2); hold on;
        else
            plot([75 100], [right_mean right_mean], 'b-', 'LineWidth', 2); hold on;
        end
        xlim([-10 110]);
        ymin = min([south_ratios north_ratios], [], 'all', 'omitmissing')*0.8; ymax = max([south_ratios north_ratios], [], 'all', 'omitmissing')*1.2;
        ymax = max(abs(log10(ymax)), abs(log10(ymin)));
        ylim([10^-ymax 10^ymax]);
        box on;
        set(gca, 'xtick', 0:25:100, 'yscale', 'log', 'ytick', 4.^(-10:10));

        if k1 == 1
            title("N cutoff=" + cutoffs(k2) + "\circ");
        end
        if k2 == 1
            ylabel("S cutoff=" + (-cutoffs(k1)) + "\circ");
        end
        if k2 > 1
            set(gca, 'yticklabel', []);
        end
        if k1 < length(cutoffs)
            set(gca, 'xticklabel', []);
        end
    end
end

end