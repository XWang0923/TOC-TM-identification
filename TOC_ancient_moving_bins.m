close all;
format longG
if ~exist('data', 'var')
    data = readtable("DM-SED TOC.xlsx");
    % Sort the imported data
    [N, M] = size(data);
    ModLat = zeros(N, 1); ModLng = zeros(N, 1); % modern coordinate
    PalLat = zeros(N, 1); PalLng = zeros(N, 1); % paleo coordinate
    for i = 1:N
        ModLat(i) = coordiCalcul(char(data{i, 2})); % modern latitude
        ModLng(i) = coordiCalcul(char(data{i, 3})); % modern longitude
        PalLat(i) = coordiCalcul(char(data{i, 4})); % paleo latitude
        PalLng(i) = coordiCalcul(char(data{i, 5})); % paleo longitude
    end
    SamAge = data{:, 1}; % interpreted age (Ma)
    TOC = data{:, 7};

    SamType = data{:, 6}; % lithology
end

% Lithology classification
% Can manually input from Excel, if required
% See details in lithology.xlsx 
% !!!!!!!!!!
% clearvars LithBins
% !!!!!!!!!!
if ~exist('LithBins', 'var')
    lithClass = readtable("lithology.xlsx", "Sheet", "All", ReadVariableNames=false);
    LithBins.dolo = lithClass{lithClass{:, 2}=="dolo", 1};
    LithBins.lime = lithClass{lithClass{:, 2}=="lime", 1};
    LithBins.carb = lithClass{lithClass{:, 2}=="carb", 1};
    LithBins.glac = lithClass{lithClass{:, 2}=="glac", 1};
    LithBins.shale = lithClass{lithClass{:, 2}=="shale", 1};
    LithBins.clay = lithClass{lithClass{:, 2}=="clay", 1};
    LithBins.sand = lithClass{lithClass{:, 2}=="sand", 1};
    LithBins.sil = lithClass{lithClass{:, 2}=="sil", 1};
    LithBins.mud = lithClass{lithClass{:, 2}=="mud", 1};
    LithBins.grain = lithClass{lithClass{:, 2}=="grain", 1};
end
% disp(["N="+num2str(length(Mn)), "Phan="+num2str(sum(SamAge<=540)), "PreC="+num2str(sum(SamAge>540))]);
clearvars -except data ModLat ModLng PalLat PalLng SamAge SamType TOC LithBins
close all;

% ageBins = [2500 1600 1000 720 635 538.8 485.4 443.8 419.2 358.9 298.9 251.902 201.4 145];
%               P_p  M_p   T   C   E     C     O     S     D     C     P       T     J
ageBins = [635 443.8 298.9 201.4];

delete results_ancient_scan_TOC.pdf
delete(gcp('nocreate'));
mypar = parpool;

for i = 1:2
    figure('Visible', 'off');
    tiledlayout(2, 1, "TileSpacing", "compact");
    idx = LithCrit(SamType, "all", LithBins) & (ageBins(i+1)<=SamAge) & (SamAge<=ageBins(i));
    Y = TOC(idx);
    X = PalLat(idx);

    [X, idx] = sort(X);
    Y = Y(idx);
    N = length(Y);

    nexttile;
    xx = floor(min(X)):ceil(max(X)); xx = xx';
    prc = zeros(length(xx), 10);
    trend_width = 5;
    for j = 1:length(xx)
        idx = abs(X-xx(j)) <= trend_width/2;
        prc(j, :) = [prctile(Y(idx), [0 5 10 25 50 75 90 95 100]) mean(Y(idx), 'omitmissing')];
    end
    for j = 1:4
        x = [xx; flipud(xx)];
        y = [prc(:, j); flipud(prc(:, 10-j))];
        x(isnan(y)) = []; y(isnan(y)) = [];
        fill(x, y, [1 1 1] * (1-0.2*j), 'EdgeColor', 'none'); hold on;
    end
    x = xx; y = prc(:, 5);
    x(isnan(y)) = []; y(isnan(y)) = [];
    plot(x, y, 'k-', 'LineWidth', 1.5); hold on;
    set(gca, 'yscale', 'log', 'xtick', -90:30:90);
    legend(["min-max", "5th-95th", "10th-90th", "25th-75th", "median"]);
    xlim([-90 90]);
    title(ageBins(i) + "-" + ageBins(i+1) + " Ma, N=" + N + " (window=" + trend_width + "\circ)");
    ylabel("TOC (wt.%)");

    nexttile;
    tar = [8 6 5];
    for j = 1:3
        x = xx;
        y = prc(:, tar(j));
        x(isnan(y)) = []; y(isnan(y)) = [];
        plot(x, y, '-', 'LineWidth', 1, 'Color', [1 1 1] * 0.4*(3-j)); hold on;
    end
    x = xx; y = prc(:, 10);
    x(isnan(y)) = []; y(isnan(y)) = [];
    plot(x, y, 'r-', 'LineWidth', 1.5); hold on;
    set(gca, 'xtick', -60:20:60);
    legend(["95th", "75th", "median", "average"]);
    xlim([-60 60]);
    xlabel("Paleolatitude (\circ)"); ylabel("TOC (wt.%)");

    exportgraphics(gcf, "results_ancient_scan_TOC.pdf",'Append',true);
    
    tar = 1:10;
    cutoffs = 15:5:35;
    yylabels = ["min", "5th", "10th", "25th", "median", "75th", "90th", "95th", "max", "average"];
    south_ratios = nan(length(cutoffs), length(cutoffs), length(tar));
    north_ratios = nan(length(cutoffs), length(cutoffs), length(tar));
    
    for j = 1:length(tar)
        disp("Working on: " + yylabels(j));
        x = xx; y = prc(:, j);
        x(isnan(y)) = []; y(isnan(y)) = [];
        [t, sratios, nratios] = regress_gauss1(x, y, cutoffs);
        south_ratios(:, :, j) = sratios;
        north_ratios(:, :, j) = nratios;
        title(t, ageBins(i) + "-" + ageBins(i+1) + " Ma, target=" + yylabels(j));
        xlabel(t, "Paleolatitude (\circ)");
        ylabel(t, "TOC (wt.%)");
        exportgraphics(gcf, "results_ancient_scan_TOC.pdf",'Append',true);
    end

    [t, ret] = regress_cross_comparison(cutoffs, south_ratios, north_ratios);
    gap = nan(length(cutoffs), 1);
    writematrix([ret(:, :, 1) gap ret(:, :, 2) gap ret(:, :, 3); gap' nan gap' nan gap'; ...
        ret(:, :, 4) gap ret(:, :, 5) gap ret(:, :, 6)], "TOC_peak_comparison.xlsx", "Sheet", "ancient_scan"+i);
    xlabel(t, "Percentiles");
    ylabel(t, "Peak height ratio");
    title(t, ageBins(i) + "-" + ageBins(i+1) + " Ma, Tropic TOC peak vs. South("+char(9661)+")/North("+char(9651)+") TOC peak");
    exportgraphics(gcf, "results_ancient_scan_TOC.pdf",'Append',true);
end

delete(mypar);
winopen("results_ancient_scan_TOC.pdf");




function ret=LithCrit(x, name, LithBins) % Lithology constraint
dolo = matches(x, LithBins.dolo);
lime = matches(x, LithBins.lime);
carb = matches(x, LithBins.carb) | dolo | lime;
shale = matches(x, LithBins.shale);
mud = matches(x, LithBins.mud);
clay = matches(x, LithBins.clay);
fineSil = shale | mud | clay;
grain = matches(x, LithBins.grain);
silici = matches(x, LithBins.sand) | fineSil | grain;
glac = matches(x, LithBins.glac);
unclass = x=="";
all = ones(length(x), 1);
ret = eval(name);
end

function ret=coordiCalcul(x) % convert string to double
if x == "0"
    ret = 0;
else
    id1 = strfind(x, "°");
    ret = str2double(x(1:id1-1));
    idx = strfind(x, "'");
    id2 = idx(1); id3 = idx(2);
    if ret >= 0
        ret = ret + str2double(x(id1+1:id2-1)) / 60;
        ret = ret + str2double(x(id2+1:id3-1)) / 60 / 60;
    else
        ret = ret - str2double(x(id1+1:id2-1)) / 60;
        ret = ret - str2double(x(id2+1:id3-1)) / 60 / 60;
    end
end
end