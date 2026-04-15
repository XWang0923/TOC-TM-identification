age1 = 12;
suffix = ["Mesoprot", "Toni", "Cryo", "Edia", "Cam", "Ord", "Sil", "Dev", "Car", "Per", "Tri"];


for j = 1:11

    age1 = j + 1;

    idx = (ageBins(age1+1)<=SamAge) & (SamAge<=ageBins(age1)) & LithCrit(SamType, "shale", LithBins);
    
    plat = PalLat(idx); plng = PalLng(idx);
    Y = Mo(idx) ./ TOC(idx);
    plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    fil = "ISEA-shale\Mo_TOC_"+suffix(j)+".xlsx";
    writetable(table(plat, plng, Y), fil);
    
    plat = PalLat(idx); plng = PalLng(idx);
    Y = MoEF(idx) ./ UEF(idx);
    plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    fil = "ISEA-shale\MoEF_UEF_"+suffix(j)+".xlsx";
    writetable(table(plat, plng, Y), fil);
    
    plat = PalLat(idx); plng = PalLng(idx);
    Y = TOC(idx) ./ P(idx) / 12 * 31;
    plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    fil = "ISEA-shale\TOC_P_"+suffix(j)+".xlsx";
    writetable(table(plat, plng, Y), fil);
    
    % plat = PalLat(idx); plng = PalLng(idx);
    % Y = TOC(idx);
    % plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    % plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    % fil = "ISEA-shale\TOC_"+suffix(j)+".xlsx";
    % writetable(table(plat, plng, Y), fil);
    % 
    % plat = PalLat(idx); plng = PalLng(idx);
    % Y = TS(idx);
    % plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    % plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    % fil = "ISEA-shale\TS_"+suffix(j)+".xlsx";
    % writetable(table(plat, plng, Y), fil);
    % 
    % plat = PalLat(idx); plng = PalLng(idx);
    % Y = P(idx);
    % plat(isnan(Y)) = []; plng(isnan(Y)) = []; Y(isnan(Y)) = [];
    % plat(Y<=0) = []; plng(Y<=0) = []; Y(Y<=0) = [];
    % fil = "ISEA-shale\P_"+suffix(j)+".xlsx";
    % writetable(table(plat, plng, Y), fil);

end

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