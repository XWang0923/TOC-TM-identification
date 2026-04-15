% Read data from Excel
% Only do once!! Tooooo much time!
% Can manually input from Excel, if required
% !!!!!!!!!!
% clear all;
% !!!!!!!!!!
format longG
if ~exist('data', 'var')
    data = readtable("Geochemistry(Oct152025).xlsx");
    % Sort the imported data
    [N, M] = size(data);
    ModLat = zeros(N, 1); ModLng = zeros(N, 1); % modern coordinate
    PalLat = zeros(N, 1); PalLng = zeros(N, 1); % paleo coordinate
    for i = 1:N
        ModLat(i) = coordiCalcul(char(data{i, 26})); % modern latitude
        ModLng(i) = coordiCalcul(char(data{i, 27})); % modern longitude
        PalLat(i) = coordiCalcul(char(data{i, 28})); % paleo latitude
        PalLng(i) = coordiCalcul(char(data{i, 29})); % paleo longitude
    end
    SamAge = data{:, 16}; % interpreted age (Ma)
    Cr = str2double(data{:, 94});
    Mo = str2double(data{:, 111});
    U = str2double(data{:, 131});
    V = str2double(data{:, 132});
    Re = str2double(data{:, 118});
    Al = data{:, 72};
    Ce = str2double(data{:, 92});
    Pr = str2double(data{:, 116});
    Nd = str2double(data{:, 113});
    TOC = data{:, 38};
    P = data{:, 79};
    TS = str2double(data{:, 80});
    TN = str2double(data{:, 77});

    % Modern sites
    emptyPal = (PalLat==0) & (PalLng==0) & (SamAge<1);
    PalLat(emptyPal) = ModLat(emptyPal);
    PalLng(emptyPal) = ModLng(emptyPal);
    % Ancient sites, without paleo-coordinates/Mn data
    emptyPal = ((PalLat==0) & (PalLng==0)) | (isnan(Cr) & isnan(Mo) & isnan(U) & isnan(V) & isnan(Al) & isnan(TOC) & isnan(P) & isnan(TS) & isnan(TN));
    data(emptyPal, :) = []; % remove from table

    emptyPal = ~emptyPal; % remove from existing data
    ModLat = ModLat(emptyPal); ModLng = ModLng(emptyPal);
    PalLat = PalLat(emptyPal); PalLng = PalLng(emptyPal);
    SamAge = SamAge(emptyPal);
    Cr = Cr(emptyPal);
    Mo = Mo(emptyPal);
    U = U(emptyPal);
    V = V(emptyPal);
    Re = Re(emptyPal);
    Al = Al(emptyPal);
    Ce = Ce(emptyPal);
    Pr = Pr(emptyPal);
    Nd = Nd(emptyPal);
    TOC = TOC(emptyPal);
    P = P(emptyPal);
    TS = TS(emptyPal);
    TN = TN(emptyPal);

    SamType = data{:, 32}; % lithology
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
clearvars -except data ModLat ModLng PalLat PalLng SamAge SamType Cr Mo U V Re Al Ce Pr Nd TOC P TS TN LithBins


close all;

% PAAS
Al0 = 18.9 * 27 / (27 + 24);
Cr0 = 110; Mo0 = 1.0; U0 = 3.1; V0 = 150;
Ce0 = 80; Pr0 = 8.9; Nd0 = 32;
% UCC
Re0_UCC = 0.198e-3; Al0_UCC = 15.40 * 27 / (27 + 24);
Cr0_UCC = 92; Mo0_UCC = 1.1; U0_UCC = 2.7; V0_UCC = 97;
Ce0_UCC = 63; Pr0_UCC = 7.1; Nd0_UCC = 27;
% enrichment factors
CeN = Ce / Ce0; PrN = Pr / Pr0; NdN = Nd / Nd0;
CrEF = Cr ./ Al ./ (Cr0 / Al0);
MoEF = Mo ./ Al ./ (Mo0 / Al0);
UEF = U ./ Al ./ (U0 / Al0);
VEF = V ./ Al ./ (V0 / Al0);
ReEF = Re ./ Al ./ (Re0_UCC / Al0_UCC);
Ce_anom = CeN ./ (PrN .* PrN ./ NdN);
bwidth = 10; % grouping - latitude range
zz = bwidth/2:bwidth:90;
zz = [-fliplr(zz) zz];

ageBins = [2500 1600 1000 720 635 538.8 485.4 443.8 419.2 358.9 298.9 251.902 201.4 145];
% ageBins = [1600 1000 540 400 250 100];
% ageBins = [2500 538.5 419.2 145];
% ageBins = [2500 1600 1400 1200 1000 720 635 575 538.8 521 509 497 485.4 470.0 458.4  443.8 433.4 427.4 423 419.2 393.3 382.7 358.9 323.2 298.9 273.01 259.51 251.902 247.2 237 201.4 174.7 161.5 145 100.5 66 23.03 2.58];
N_bins = length(ageBins) - 1;

figure;

%%
grad = zeros(N_bins, 7);
R2 = zeros(N_bins, 1);
sample_cnt = zeros(N_bins, 1);
bin_cnt = zeros(N_bins, 1);
bin_rec = zeros(N_bins, length(zz));
tropic_grad = zeros(N_bins, 1);
tropic_R2 = zeros(N_bins, 1);
for i = 1:N_bins
    nexttile;

    idx = LithCrit(SamType, "all", LithBins) & (ageBins(i+1)<=SamAge) & (SamAge<=ageBins(i));
    Y = TOC(idx) ./ TN(idx) / 12 * 14;
    % Y(Y<0) = NaN;
    plat = PalLat(idx);
    plat(isnan(Y)) = []; Y(isnan(Y)) = [];
    plat(~Y) = []; Y(~Y) = [];
    g = (round(plat/bwidth+0.5)-0.5)*bwidth;
    % g = abs(g); % Latitude-flipping
    

    boxg = g; boxy = Y;
    plt = boxchart(boxg, boxy, 'MarkerStyle', 'x', 'BoxFaceAlpha', 0.6); hold on;
    plt.BoxWidth = 5;
    plt.JitterOutliers = "on";
    xlim([-90 90]);
    box on;
    set(gca, 'yscale', 'log', 'xtick', -90:30:90);
    title(ageBins(i) + "-" + ageBins(i+1) + "Ma, N=" + sum(~isnan(Y)));
    sample_cnt(i) = sum(~isnan(Y));

    % Bootstrap
    xx = unique(g);
    bin_cnt(i) = length(xx);
    iter = 101;
    res = zeros(iter, 4);
    nn = length(Y);
    for k = 1:iter
        % pick = randi(sample_cnt(i), sample_cnt(i), 1);
        pick = randperm(sample_cnt(i), max(50, round(0.05*sample_cnt(i))));
        % prob = 3 - floor(abs(g)/30);
        % prob = prob .^ 0.5;
        % pick = datasample(1:nn, rr, 'Weights', prob, 'Replace', false);
        Y_sub = Y(pick); g_sub = g(pick);
        yy = 0 * xx;
        for j = 1:bin_cnt(i)
            yy(j) = prctile(Y_sub(g_sub==xx(j)), 95);
            if ~isnan(yy(j))
                bin_rec(i, zz==xx(j)) = bin_rec(i, zz==xx(j)) + 1;
            end
        end
        xx_sub = xx(~isnan(yy)); yy_sub = yy(~isnan(yy));
        plot(xx_sub, yy_sub, 'Color', [[1 1 1]*0.4 0.1]);
        yy_log = log10(yy_sub);
        [b, bint, r, rint, stats] = regress(yy_log, [ones(length(xx_sub), 1) xx_sub xx_sub.*xx_sub]);
        res(k, 1:3) = b'; res(k, 4) = stats(1);
    end

    b = mean(res, 1, 'omitmissing');
    b_sd = std(res, 1, 'omitmissing');
    grad(i, :) = [prctile(res(:, 3), [10 25 50 70 90]) b(3) b_sd(3)]; R2(i) = b(4);

    idx = bin_rec(i, :) >= iter*0.6;
    xx = min(zz(idx)):bwidth:max(zz(idx));
    yy = b(1) + b(2)*xx + b(3).*xx.*xx;
    yy = 10 .^ yy;
    plot(xx, yy, 'r', 'LineWidth', 3); hold on;
    % [xx' yy']

    xx = unique(g);
    xx(xx<-40 | xx>40) = [];
    yy = 0 * xx;
    for j = 1:length(xx)
        yy(j) = prctile(Y(g==xx(j)), 95);
    end
    yy = log10(yy);
    [b, ~, ~, ~, stats] = regress(yy, [ones(length(xx), 1) xx xx.*xx]);
    tropic_grad(i) = b(3);
    tropic_R2(i) = stats(1);

    xx = min(xx):bwidth:max(xx);
    yy = b(1) + b(2) * xx + b(3) .* xx .* xx;
    yy = 10 .^ yy;
    plot(xx, yy, 'k', 'LineWidth', 2); hold on;

    % break;
end

% [xx' yy']

xx = -85:10:85;
output = nan(length(boxy), length(xx));
for i = 1:length(xx)
    yy = Y(g==xx(i));
    nn = length(yy);
    output(1:nn+1, i) = [xx(i); yy];
end
writematrix(output, "boxchart.xlsx");

% figure;
% for i = 1:N_bins
%     nexttile;
%     idx = ~~bin_rec(i, :);
%     scatter(zz(idx), bin_rec(i, idx), 'b', 'filled'); hold on;
%     plot(zz, 0*zz+iter*0.6, 'r-'); hold on;
%     box on;
%     ylim([-10 iter+10]);
%     title(ageBins(i) + "-" + ageBins(i+1) + "Ma");
% end

samLimit = 250;
binLimit = 8;

figure;
X = (ageBins(1:end-1) + ageBins(2:end)) / 2;
X = X(2:end-1);
grad = grad(2:end-1, :);
R2 = R2(2:end-1);
tropic_grad = tropic_grad(2:end-1);
tropic_R2 = tropic_R2(2:end-1);
sample_cnt = sample_cnt(2:end-1);
bin_cnt = bin_cnt(2:end-1);
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
fill([X fliplr(X)], [grad(:, 2)' fliplr(grad(:, 4)')], [1 1 1] * 0.8, 'EdgeColor', 'none'); hold on;
plot(X, 0*X, 'k--'); hold on;
plot(X, grad(:, 3), 'k-', 'LineWidth', 2); hold on;
% scatter(X, grad(:, 6)); hold on;
sameGrad = grad(:, 3) .* tropic_grad > 0;
validBin = ~(sample_cnt<samLimit | bin_cnt<binLimit);
scatter(X(sameGrad&validBin), grad(sameGrad&validBin, 3), 120, 'ro', 'filled', 'MarkerEdgeColor', 'k');
scatter(X(sameGrad&~validBin), grad(sameGrad&~validBin, 3), 120, 'bs', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
scatter(X(~sameGrad&validBin), grad(~sameGrad&validBin, 3), 120, 'bs', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
scatter(X(~sameGrad&~validBin), grad(~sameGrad&~validBin, 3), 120, 'bs', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
writematrix([X' grad(:, 3) sameGrad validBin grad(:, 2) grad(:, 4)], "output.xlsx");


ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 200:100:1000);



nexttile;
yyaxis left;
plot(X, tropic_grad, '-o'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
plot(X, 0*X, 'k--');
ylabel("Tropic gradient");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 100:100:1000);
yyaxis right;
plot(X, grad(:, 3), '-^'); hold on;
ylimit = ylim;
ylim([-1 1] * max(abs(ylimit)));
ylabel("Full Gradient");


nexttile;
plot(X, sample_cnt); hold on;
plot(X, samLimit+0*X, '--'); hold on;
ylabel("Sample count");
xlim([200 1500]);
set(gca, 'xscale', 'log', 'xtick', 100:100:1000);

nexttile;
plot(X, bin_cnt); hold on;
plot(X, binLimit-0.5+0*X, '--');
ylabel("Latitudinal bin count");
xlim([200 1500]); ylim([5 14]);
set(gca, 'xscale', 'log', 'xtick', 100:100:1000);


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
