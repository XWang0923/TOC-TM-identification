suffix = ["Mesoprot", "Toni", "Cryo", "Edia", "Cam", "Ord", "Sil", "Dev", "Car", "Per", "Tri"];
ageBins = [2500 1600 1000 720 635 538.8 485.4 443.8 419.2 358.9 298.9 251.902 201.4 145];

close all;


res = zeros(11, 2);
sam_cnt = zeros(11, 1);
figure;
for j = 6:11

    age1 = j + 1;

    idx = (ageBins(age1+1)<=SamAge) & (SamAge<=ageBins(age1)) & LithCrit(SamType, "all", LithBins);
    Y = VEF(idx);
    Y(Y<1) = NaN;

    plat = PalLat(idx);
    plat(isnan(Y)) = []; Y(isnan(Y)) = [];
    plat(Y<=0) = []; Y(Y<=0) = [];

    sam_cnt(j) = sum(~isnan(Y));
    nn = length(plat);
    [xx, I] = sort(plat);
    yy = Y(I);
    
    N = 30;
    mm = floor(nn/N);
    M = nn - N * mm;
    front = 1;
    x = zeros(N, 1); y = zeros(N, 1);
    for i = 1:N
        if i <= M
            tail = front + mm;
        else
            tail = front + mm - 1;
        end
        t1 = xx(front:tail);
        t2 = yy(front:tail);
        x(i) = mean(t1); y(i) = prctile(t2, 95);
        front = tail + 1;
    end
    

    nexttile;
    scatter(x, y, 'LineWidth', 1); hold on;
    [x, y]

    y = log10(y);
    b = regress(y, [ones(N, 1) x x.*x]);
    res(j, 1) = b(3);
    xx = floor(min(x)):ceil(max(x));
    yy = b(1) + b(2) * xx + b(3) * xx .* xx;
    yy = 10 .^ yy;
    plot(xx, yy, 'r-', 'LineWidth', 3); hold on;
    set(gca, 'yscale', 'log');
    xlim([-90 90]);
    [xx', yy']

    y = y(abs(x)<=40); x = x(abs(x)<=40);
    b = regress(y, [ones(length(x), 1) x x.*x]);
    res(j, 2) = b(3);
    xx = floor(min(x)):ceil(max(x));
    yy = b(1) + b(2) * xx + b(3) * xx .* xx;
    yy = 10 .^ yy;
    plot(xx, yy, 'k-', 'LineWidth', 2); hold on;
    [xx', yy']

    title(suffix(j)); box on;
    break;
    
end

ret = nan(11, 3);
cc = sam_cnt>250;
for i = 1:11
    if res(i, 1)*res(i, 2)>0 && cc(i)
        ret(i, :) = [res(i, 1) res(i, 1) NaN];
    else
        ret(i, :) = [res(i, 1) NaN res(i, 1)];
    end
end

% ret

X = (ageBins(1:end-1) + ageBins(2:end)) / 2;
X = X(2:end-1);

figure;
nexttile;
plot(X, res(:, 1), "k-"); hold on;
scatter(X, ret(:, 2), 120, 'ro', 'filled', 'MarkerEdgeColor', 'k'); hold on;
scatter(X, ret(:, 3), 120, 'bs', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7); hold on;
plot(X, 0*X, 'k--'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);


nexttile;
yyaxis left;
plot(X, res(:, 2), '-o'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
plot(X, 0*X, 'k--');
ylabel("Tropic gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 100:100:1000);
yyaxis right;
plot(X, res(:, 1), '-^'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Full Gradient");

nexttile;
plot(X, sam_cnt); hold on;
plot(X, 250+0*X, '--'); hold on;
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);
ylabel("Sample count");

% nexttile;


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