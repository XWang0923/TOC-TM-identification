function [t, south_ratios, north_ratios] = regress_gauss1(lat, data, cutoffs)

RGB = orderedcolors("gem");
H = rgb2hex(RGB);
H = H([1 2 4]);

gauss1 = @(p, x) p(1)*exp(-(x-p(2)).^2/(2*p(3)^2));
south_ratios = nan(length(cutoffs));
north_ratios = nan(length(cutoffs));

figure('Visible', 'off');
t = tiledlayout(length(cutoffs), length(cutoffs), "TileSpacing", "compact");
for k1 = 1:length(cutoffs)
    for k2 = 1:length(cutoffs)
        nexttile;
        x = lat; y = data;
        x(isnan(y)) = []; y(isnan(y)) = [];
        y = y(abs(x)<=60); x = x(abs(x)<=60);
        plot(x, y, '-', 'LineWidth', 1, 'Color', 0.4 * [1 1 1]); hold on;
        
        fitted_ranges = nan(3, 2);
        fitted_params = nan(3, 3);
        fitted_judge = zeros(3, 1);
        fitted_peaks = nan(3, 1);
        
        boundaries = [-60 -cutoffs(k1) cutoffs(k2) 60];
        for i = 1:3
            x = lat(boundaries(i)<=lat & lat<=boundaries(i+1));
            y = data(boundaries(i)<=lat & lat<=boundaries(i+1));
            if ~isempty(x)
                fitted_ranges(i, :) = [floor(min(x)) ceil(max(x))];
                initial_guess= [10, mean(boundaries([i i+1])), 10];
                options = optimset('Display','off');
                n = length(x);
                bootstrap_n = 101;
                if bootstrap_n > 1
                    bootstrap_params = zeros(bootstrap_n, 3);
                    parfor j = 1:bootstrap_n
                        pick = randi(n, n, 1);
                        bootstrap_params(j, :) = lsqcurvefit(gauss1, initial_guess, x(pick), y(pick), [], [], options);
                    end
                    fitted_params(i, 2) = median(bootstrap_params(:, 2));
                    idx = find(fitted_params(i, 2)==bootstrap_params(:, 2), 1);
                    fitted_params(i, [1 3]) = bootstrap_params(idx, [1 3]);
                else
                    fitted_params(i, :) = lsqcurvefit(gauss1, initial_guess, x, y, [], [], options);
                end
                if boundaries(i)<=fitted_params(i, 2) && fitted_params(i, 2)<=boundaries(i+1)
                    fitted_judge(i) = 1;
                    fitted_peaks(i) = fitted_params(i, 1);
                end
            end
        end
        
        for k3 = 1:3
            x_fit = fitted_ranges(k3, 1):fitted_ranges(k3, 2);
            y_fit = gauss1(fitted_params(k3, :), x_fit);
            if fitted_judge(k3)
                plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', H(k3)); hold on;
            else
                plot(x_fit, y_fit, '--', 'LineWidth', 1.5, 'Color', H(k3)); hold on;
            end
            if k1 == 1 && k2 == 1
                [x_fit; y_fit]'
            end
        end

        south_ratios(k1, k2) = fitted_peaks(2) / fitted_peaks(1);
        north_ratios(k1, k2) = fitted_peaks(2) / fitted_peaks(3);
        xlim([-60 60]);
        set(gca, 'xtick', -60:30:60);
        
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