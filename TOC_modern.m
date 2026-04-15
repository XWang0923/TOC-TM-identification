close all; clear;

lat = ncread("toc_concentration_orig.nc", "latitude");
lng = ncread("toc_concentration_orig.nc", "longitude");
TOC = ncread("toc_concentration_orig.nc", "toc_concentration");
TOC(~TOC) = NaN;
TOC = TOC';

delete results_modern_TOC.pdf
delete(gcp('nocreate'));
mypar = parpool;

figure('Visible', 'off');
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
h = imagesc(lng, lat, TOC); hold on;
set(gca, 'ydir', 'normal');
set(h, 'alphadata', ~isnan(TOC));
clim([0 5]);
cb = colorbar;
ylabel(cb, "TOC (wt.%)");
xlabel("Longitude (\circ)"); ylabel("Latitude (\circ)");
title("Global TOC map");

TOC(lat>=30&lat<=66, lng>=10&lng<=60) = NaN;
TOC(lat<-60|lat>60, :) = NaN;

nexttile;
h = imagesc(lng, lat, TOC); hold on;
set(gca, 'ydir', 'normal');
set(h, 'alphadata', ~isnan(TOC));
clim([0 5]);
colorbar;
cb = colorbar;
ylabel(cb, "TOC (wt.%)");
xlabel("Longitude (\circ)"); ylabel("Latitude (\circ)");
title("Global TOC map (without large restricted basins/polar regions)");

nexttile;
res0 = prctile(TOC, 0, 2);
res5 = prctile(TOC, 5, 2);
res10 = prctile(TOC, 10, 2);
res25 = prctile(TOC, 25, 2);
res50 = prctile(TOC, 50, 2);
res75 = prctile(TOC, 75, 2);
res90 = prctile(TOC, 90, 2);
res95 = prctile(TOC, 95, 2);
res100 = prctile(TOC, 100, 2);

lat(isnan(res0)) = [];
res0(isnan(res0)) = [];
res5(isnan(res5)) = [];
res10(isnan(res10)) = [];
res25(isnan(res25)) = [];
res50(isnan(res50)) = [];
res75(isnan(res75)) = [];
res90(isnan(res90)) = [];
res95(isnan(res95)) = [];
res100(isnan(res100)) = [];

fill([lat; flipud(lat)], [res0; flipud(res100)], [1 1 1] * 0.8, 'EdgeColor', 'none'); hold on;
fill([lat; flipud(lat)], [res5; flipud(res95)], [1 1 1] * 0.6, 'EdgeColor', 'none'); hold on;
fill([lat; flipud(lat)], [res10; flipud(res90)], [1 1 1] * 0.4, 'EdgeColor', 'none'); hold on;
fill([lat; flipud(lat)], [res25; flipud(res75)], [1 1 1] * 0.2, 'EdgeColor', 'none'); hold on;
plot(lat, res50, 'k-', 'LineWidth', 1.5); hold on;
set(gca, 'yscale', 'log');
legend(["min-max", "5th-95th", "10th-90th", "25th-75th", "median"]);
xlabel("Latitude (\circ)"); ylabel("TOC (wt.%)");
title("Latitudinal TOC distribution");

nexttile;
res = mean(TOC, 2, "omitmissing");
res(isnan(res)) = [];
plot(lat, res95, '-', 'LineWidth', 1, 'Color', [1 1 1] * 0.8); hold on;
plot(lat, res75, '-', 'LineWidth', 1, 'Color', [1 1 1] * 0.4); hold on;
plot(lat, res50, 'k-', 'LineWidth', 1); hold on;
plot(lat, res, 'r-', 'LineWidth', 1.5); hold on;
legend(["95th", "75th", "median", "average"]);
xlabel("Latitude (\circ)"); ylabel("TOC (wt.%)");
title("Average vs. percentiles");

exportgraphics(gcf, "results_modern_TOC.pdf",'Append',true);


gauss1 = @(p, x) p(1)*exp(-(x-p(2)).^2/(2*p(3)^2));

prc = [res0 res5 res10 res25 res50 res75 res90 res95 res100 res];
yylabels = ["min", "5th", "10th", "25th", "median", "75th", "90th", "95th", "max", "average"];


tar = 1:10;
cutoffs = 15:5:35;
south_ratios = nan(length(cutoffs), length(cutoffs), 10);
north_ratios = nan(length(cutoffs), length(cutoffs), 10);
for j = 1:length(tar)
    disp("Working on: " + yylabels(j));
    [t, sratios, nratios] = regress_gauss1(lat, prc(:, j), cutoffs);
    south_ratios(:, :, j) = sratios;
    north_ratios(:, :, j) = nratios;
    title(t, yylabels(tar(j)) + " in modern latitudinal TOC distribution");
    xlabel(t, "Paleolatitude (\circ)");
    ylabel(t, "TOC (wt.%)");
    exportgraphics(gcf, "results_modern_TOC.pdf",'Append',true);
end

[t, ret] = regress_cross_comparison(cutoffs, south_ratios, north_ratios);
gap = nan(length(cutoffs), 1);
writematrix([ret(:, :, 1) gap ret(:, :, 2) gap ret(:, :, 3); gap' nan gap' nan gap'; ...
    ret(:, :, 4) gap ret(:, :, 5) gap ret(:, :, 6)], "TOC_peak_comparison.xlsx", "Sheet", "modern");
xlabel(t, "Percentiles");
ylabel(t, "Peak height ratio");
title(t, "Modern, Tropic TOC peak vs. South("+char(9661)+")/North("+char(9651)+") TOC peak");
exportgraphics(gcf, "results_modern_TOC.pdf",'Append',true);

delete(mypar);
winopen("results_modern_TOC.pdf");
