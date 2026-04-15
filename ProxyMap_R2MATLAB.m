suffix = ["Mesoprot", "Toni", "Cryo", "Edia", "Cam", "Ord", "Sil", "Dev", "Car", "Per", "Tri"];
ageBins = [2500 1600 1000 720 635 538.8 485.4 443.8 419.2 358.9 298.9 251.902 201.4 145];

close all;

ret1 = zeros(11, 1);
ret2 = zeros(11, 1);
sam_cnt = zeros(11, 1);
grid_cnt = zeros(11, 1);
figure;
for j = 6:11

    fil = "ISEA-all\VEF_" + suffix(j) + ".xlsx";
    
    
    data1 = readtable(fil, "Sheet", "Sheet2");
    data2 = readtable(fil, "Sheet", "Sheet3");
    data3 = readtable(fil, "Sheet", "Sheet4");
    
    plat = data1{:, 2};
    plng = data1{:, 3};
    Y = data1{:, 4};
    loc = data1{:, 5};
    
    Y(Y<1) = NaN;

    plat(isnan(Y)) = []; plng(isnan(Y)) = [];
    loc(isnan(Y)) = []; Y(isnan(Y)) = [];

    sam_cnt(j) = length(Y);

    poly_v = data2{:, 1};
    polyID = data2{:, 2};
    
    grid_cnt(j) = length(polyID);
    
    % nexttile;
    geo_v = data3{:, 1};
    geo_cor = data3{:, 2};
    N = length(geo_cor);
    poly_lng = zeros(N, 1);
    poly_lat = zeros(N, 1);
    for i = 1:length(geo_v)
        str = string(geo_cor(i));
        num_str = regexp(str,'\-?\d*\.?\d*','match');
        nn = length(num_str) / 2;
        xx = str2double(num_str(1:nn));
        yy = str2double(num_str(nn+1:end));
        % plot(xx, yy, 'k-'); hold on;
        poly_lng(i) = mean(xx(1:end-1)); poly_lat(i) = mean(yy(1:end-1));
    end
    
    tot_cnt = zeros(N, 1);
    res = zeros(N, 4);
    for i = 1:N
        idx = polyID(i)==loc;
        tot_cnt(i) = sum(idx);
        yy = Y(idx);
        res(i, :) = [mean(yy) prctile(yy, [50 75 95])];
    end
    
    % scatter(poly_lng, poly_lat, 15, tot_cnt, 'filled'); hold on;
    % colorbar;
    
    nexttile;
    xx = poly_lat; yy = res(:, 4);
    scatter(xx, yy, 'LineWidth', 1); hold on;
    [xx, yy]
    yy = log10(yy);
    [b,~,~,~,~] = regress(yy, [ones(length(xx), 1) xx xx.*xx]);
    xr = round(min(xx)):round(max(xx));
    yr = b(1) + b(2) * xr + b(3) .* xr .* xr;
    yr = 10 .^ yr;
    plot(xr, yr, 'r-', 'linewidth', 3); hold on;
    set(gca, 'yscale', 'log');
    [xr', yr']
    
    ret1(j) = b(3);
    
    xx = poly_lat; yy = res(:, 4);
    yy = yy(abs(xx)<=40); xx = xx(abs(xx)<=40);
    yy = log10(yy);
    [b,~,~,~,~] = regress(yy, [ones(length(xx), 1) xx xx.*xx]);
    xr = round(min(xx)):round(max(xx));
    yr = b(1) + b(2) * xr + b(3) .* xr .* xr;
    yr = 10 .^ yr;
    plot(xr, yr, 'k', 'LineWidth', 2); hold on;
    [xr', yr']

    box on;
    xlim([-90 90]);
    title(suffix(j));

    ret2(j) = b(3);

    break;
end

X = (ageBins(1:end-1) + ageBins(2:end)) / 2;
X = X(2:end-1);

figure;
tiledlayout(2, 2, 'TileSpacing', 'compact');



ret = nan(11, 3);
ret(:, 1) = ret1;
id1 = ret1.*ret2>0;
id2 = sam_cnt>250;
id3 = grid_cnt>20;
id = id1 & id2 & id3;
for i = 1:11
    if id1(i) && id2(i) && id3(i)
        ret(i, 2) = ret1(i);
    else
        ret(i, 3) = ret1(i);
    end
end

% ret

nexttile;
plot(X, ret1, "k-"); hold on;
scatter(X(id), ret1(id), 120, 'ro', 'filled', 'MarkerEdgeColor', 'k'); hold on;
scatter(X(~id), ret1(~id), 120, 'bs', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7); hold on;
plot(X, 0*X, 'k--'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);


nexttile;
yyaxis left;
plot(X, ret2, '-o'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
plot(X, 0*X, 'k--');
ylabel("Tropic gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 100:100:1000);
yyaxis right;
plot(X, ret1, '-^'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Full Gradient");

nexttile;
plot(X, sam_cnt); hold on;
plot(X, 250+0*X, '--'); hold on;
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);
ylabel("Sample count");

nexttile;
plot(X, grid_cnt); hold on;
plot(X, 20+0*X, '--'); hold on;
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);
ylabel("Cell count");