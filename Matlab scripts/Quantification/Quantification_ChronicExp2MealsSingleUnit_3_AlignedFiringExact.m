%% Housekeeping
close all;clearvars;clc;
%%
dir_info = dir("*_recordingInfo.mat");
N_recordings = length(dir_info);
load('../Analysis_ENS_ChronicExp2MealTimes.mat');
%%
n_cells = zeros(1, N_recordings);

% binSize_vec = [1, 2, 5, 10, 20];
binSize_vec = [1, 3, 5, 15];
n_binSize = length(binSize_vec);

RatID = [];
ExpDay = [];

histAll = [];
histMealAll = [];
histBeforeAll = [];
histDuringAll = [];
histAfterAll = [];
zHistAll = [];
zHistMealAll = [];
zHistBeforeAll = [];
zHistDuringAll = [];
zHistAfterAll = [];
mzBeforeAll = [];
mzDuringAll = [];
mzAfterAll = [];
lmBeforeAll = [];
lmDuringAll = [];
lmAfterAll = [];
slpBeforeAll = [];
slpDuringAll = [];
slpAfterAll = [];
fitBeforeAll = [];
fitDuringAll = [];
fitAfterAll = [];

for idxFile = 1:N_recordings

    load(dir_info(idxFile).name);
    nameFile = recordingInfo.recordingName;
    ratIDFile = str2double(nameFile(end));
    expDayFile = str2double(nameFile(24:end-5));

    MealTimeMin = tblMealTime.MealTimeMinutes{ strcmp(nameFile, ...
        tblMealTime.RecordingName)};
    MealTimeExactSec = tblMealTime.MealTimeExactSec{ strcmp(nameFile, ...
        tblMealTime.RecordingName)};
    FirstMealStart = MealTimeExactSec(1);
    FirstMealStop = FirstMealStart + 60;

    load([nameFile, '_ManualSingleUnit.mat'], 'su_t');
    nCell = length(su_t);
    n_cells(idxFile) = nCell;
    
    histFile = cell(nCell, n_binSize);
    histMealFile = cell(nCell, n_binSize);
    histBeforeFile = cell(nCell, n_binSize);
    histDuringFile = cell(nCell, n_binSize);
    histAfterFile = cell(nCell, n_binSize);
    zHistFile = cell(nCell, n_binSize);
    zHistMealFile = cell(nCell, n_binSize);
    zHistBeforeFile = cell(nCell, n_binSize);
    zHistDuringFile = cell(nCell, n_binSize);
    zHistAfterFile = cell(nCell, n_binSize);
    mzBeforeFile = zeros(nCell, n_binSize);
    mzDuringFile = zeros(nCell, n_binSize);
    mzAfterFile = zeros(nCell, n_binSize);
    lmBeforeFile = cell(nCell, n_binSize);
    lmDuringFile = cell(nCell, n_binSize);
    lmAfterFile = cell(nCell, n_binSize);
    slpBeforeFile = zeros(nCell, n_binSize);
    slpDuringFile = zeros(nCell, n_binSize);
    slpAfterFile = zeros(nCell, n_binSize);
    fitBeforeFile = cell(nCell, n_binSize);
    fitDuringFile = cell(nCell, n_binSize);
    fitAfterFile = cell(nCell, n_binSize);

    tMealVec = cell(n_binSize, 1);
    tBeforeVec = cell(n_binSize, 1);
    tDuringVec = cell(n_binSize, 1);
    tAfterVec = cell(n_binSize, 1);
    for idxSize = 1:n_binSize
        binEdges = 0:binSize_vec(idxSize):30*60;
        binCenters = (binEdges(1:end-1) + binEdges(2:end))./2;

        binsMeal = binCenters >= (FirstMealStart - 8*60) & ...
            binCenters <= (FirstMealStart + 8*60);
        tMeal = linspace(-8, 8, sum(binsMeal));
        tMealVec{idxSize} = tMeal;
        
        binsBefore = binCenters <= FirstMealStart & ...
            binCenters >= (FirstMealStart - 1*60);
        tBefore = linspace(-1, 0, sum(binsBefore));
        tBeforeVec{idxSize} = tBefore;

        binsDuring = binCenters <= FirstMealStop & ...
            binCenters >= FirstMealStart;
        tDuring = linspace(0, 1, sum(binsDuring));
        tDuringVec{idxSize} = tDuring;

        binsAfter = binCenters >= FirstMealStop & ...
            binCenters <= (FirstMealStop + 1*60);
        tAfter = linspace(0, 1, sum(binsAfter));
        tAfterVec{idxSize} = tAfter;

        for idxCell = 1:nCell

            histFile{idxCell, idxSize} = histcounts(su_t{idxCell}, ...
                binEdges) ./ binSize_vec(idxSize);
            zHistFile{idxCell, idxSize} = zscore(histFile{idxCell, idxSize});

            histMealFile{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsMeal);
            zHistMealFile{idxCell, idxSize} = zHistFile{idxCell, idxSize}( ...
                binsMeal);

            histBeforeFile{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsBefore);
            zHistBeforeFile{idxCell, idxSize} = zHistFile{idxCell, idxSize}( ...
                binsBefore);
            mzBeforeFile(idxCell, idxSize) = mean( ... 
                zHistBeforeFile{idxCell, idxSize});
            lmBeforeFile{idxCell, idxSize} = fitlm(tBefore, ...
                zHistBeforeFile{idxCell, idxSize}, 'linear');
            fitBeforeFile{idxCell, idxSize} = ...
                lmBeforeFile{idxCell, idxSize}.Fitted;
            slpBeforeFile(idxCell, idxSize) = ...
                lmBeforeFile{idxCell, idxSize}.Coefficients.Estimate(2);

            histDuringFile{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsDuring);
            zHistDuringFile{idxCell, idxSize} = zHistFile{idxCell, idxSize}( ...
                binsDuring);
            mzDuringFile(idxCell, idxSize) = mean( ... 
                zHistDuringFile{idxCell, idxSize});
            lmDuringFile{idxCell, idxSize} = fitlm(tDuring, ...
                zHistDuringFile{idxCell, idxSize}, 'linear');
            fitDuringFile{idxCell, idxSize} = ...
                lmDuringFile{idxCell, idxSize}.Fitted;
            slpDuringFile(idxCell, idxSize) = ...
                lmDuringFile{idxCell, idxSize}.Coefficients.Estimate(2);
            
            histAfterFile{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsAfter);
            zHistAfterFile{idxCell, idxSize} = zHistFile{idxCell, idxSize}( ...
                binsAfter);
            mzAfterFile(idxCell, idxSize) = mean( ...
                zHistAfterFile{idxCell, idxSize});
            lmAfterFile{idxCell, idxSize} = fitlm(tAfter, ...
                zHistAfterFile{idxCell, idxSize}, 'linear');
            fitAfterFile{idxCell, idxSize} = ...
                lmAfterFile{idxCell, idxSize}.Fitted;
            slpAfterFile(idxCell, idxSize) = ...
                lmAfterFile{idxCell, idxSize}.Coefficients.Estimate(2);

            dummy = 1;
        end
    end
    
    RatID = [RatID; ratIDFile*ones(nCell,1)];
    ExpDay = [ExpDay; expDayFile*ones(nCell,1)];
    histAll = [histAll; histFile];
    histMealAll = [histMealAll; histMealFile];
    histBeforeAll = [histBeforeAll; histBeforeFile];
    histDuringAll = [histDuringAll; histDuringFile];
    histAfterAll = [histAfterAll; histAfterFile];
    zHistAll = [zHistAll; zHistFile];
    zHistMealAll = [zHistMealAll; zHistMealFile];
    zHistBeforeAll = [zHistBeforeAll; zHistBeforeFile];
    zHistDuringAll = [zHistDuringAll; zHistDuringFile];
    zHistAfterAll = [zHistAfterAll; zHistAfterFile];
    mzBeforeAll = [mzBeforeAll; mzBeforeFile];
    mzDuringAll = [mzDuringAll; mzDuringFile];
    mzAfterAll = [mzAfterAll; mzAfterFile];
    lmBeforeAll = [lmBeforeAll; lmBeforeFile];
    lmDuringAll = [lmDuringAll; lmDuringFile];
    lmAfterAll = [lmAfterAll; lmAfterFile];
    slpBeforeAll = [slpBeforeAll; slpBeforeFile];
    slpDuringAll = [slpDuringAll; slpDuringFile];
    slpAfterAll = [slpAfterAll; slpAfterFile];
    fitBeforeAll = [fitBeforeAll; fitBeforeFile];
    fitDuringAll = [fitDuringAll; fitDuringFile];
    fitAfterAll = [fitAfterAll; fitAfterFile];

    dummy = 1;
end
N_Cells = sum(n_cells);
CellID = (1:N_Cells)';

T = table(CellID, RatID, ExpDay);
%%
for idxSize = 1:n_binSize
    T.(['MeanFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        mzBeforeAll(:, idxSize);
    T.(['MeanFiringDuring_', num2str(binSize_vec(idxSize))]) = ...
        mzDuringAll(:, idxSize);
    T.(['MeanFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        mzAfterAll(:, idxSize);
    T.(['MeanFiringDiff_', num2str(binSize_vec(idxSize))]) = ...
        mzAfterAll(:, idxSize) - mzBeforeAll(:, idxSize);

    T.(['SlopeFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        slpBeforeAll(:, idxSize);
    T.(['SlopeFiringDuring_', num2str(binSize_vec(idxSize))]) = ...
        slpDuringAll(:, idxSize);
    T.(['SlopeFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        slpAfterAll(:, idxSize);
    T.(['FittedFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        fitBeforeAll(:, idxSize);
    T.(['FittedFiringDuring_', num2str(binSize_vec(idxSize))]) = ...
        fitDuringAll(:, idxSize);
    T.(['FittedFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        fitAfterAll(:, idxSize);
    T.(['LinearModelBefore_', num2str(binSize_vec(idxSize))]) = ...
        lmBeforeAll(:, idxSize);
    T.(['LinearModelDuring_', num2str(binSize_vec(idxSize))]) = ...
        lmDuringAll(:, idxSize);
    T.(['LinearModelAfter_', num2str(binSize_vec(idxSize))]) = ...
        lmAfterAll(:, idxSize);
    T.(['FiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        histBeforeAll(:, idxSize);
    T.(['zFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        zHistBeforeAll(:, idxSize);
    T.(['FiringDuring_', num2str(binSize_vec(idxSize))]) = ...
        histDuringAll(:, idxSize);
    T.(['zFiringDuring_', num2str(binSize_vec(idxSize))]) = ...
        zHistDuringAll(:, idxSize);
    T.(['FiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        histAfterAll(:, idxSize);
    T.(['zFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        zHistAfterAll(:, idxSize);
    T.(['FiringMeal_', num2str(binSize_vec(idxSize))]) = ...
        histMealAll(:, idxSize);
    T.(['zFiringMeal_', num2str(binSize_vec(idxSize))]) = ...
        zHistMealAll(:, idxSize);
    T.(['Firing_', num2str(binSize_vec(idxSize))]) = ...
        histAll(:, idxSize);
    T.(['zFiring_', num2str(binSize_vec(idxSize))]) = ...
        zHistAll(:, idxSize);
end

%% Remove day8_rat3 and day8_rat4
% idxRemoval = (T.ExpDay == 8) & (T.RatID == 3 | T.RatID == 4);
% t = T(~idxRemoval, :);

t = T;
t.RatID = categorical(t.RatID);
t.ExpDay = categorical(t.ExpDay);
t.CellID = categorical(t.CellID);
N_Cells = height(t);
%%
fig1 = figure('Name','FiringAourndMeal', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    
    tmp_MeanFiringBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    [~, tmp_order] = sort(tmp_MeanFiringBefore, 'ascend');
    tmp_FiringMeal = vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:});
    tmp_FiringMealSorted = tmp_FiringMeal(tmp_order,:);

    subplot(1, n_binSize, i);
    imagesc(tMealVec{i}, 1:N_Cells, tmp_FiringMealSorted);hold on;grid on;box on;
    % set(gca, 'Ydir','normal');
    xlim([-4, 6]);xlabel('Time');
    ylim([1, N_Cells]);ylabel('Cells');
    title(['Ordered Firing BinSize_',num2str(binSize_vec(i)), 's'], ...
        'Interpreter','none');
    colormap jet;colorbar('southoutside');
end
%%
fig2 = figure('Name','FiringAourndMealMean', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    
    tmp_z = mean(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1);
    tmp_sem = std(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1)./sqrt(N_Cells);

    subplot(1, n_binSize, i);hold on;box on;
    shadedErrorBar(tMealVec{i}, tmp_z, tmp_sem);
    xregion(-1, 0, 'FaceColor','y');
    xregion(0, 1, 'FaceColor', 'c');
    xregion(1, 2, 'FaceColor','g');
    xlim([-4, 6]);xlabel('Time');
    title(['Mean Firing BinSize_',num2str(binSize_vec(i)), 's'], ...
        'Interpreter','none');
end
%%
KWrateP = zeros(1, n_binSize);
KWrateStat = cell(1, n_binSize);
fig3 = figure('Name','FiringAroundMealBox', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    tmp_MeanBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    tmp_MeanDuring = t.(['MeanFiringDuring_', num2str(binSize_vec(i))]);
    tmp_MeanAfter = t.(['MeanFiringAfter_', num2str(binSize_vec(i))]);
    
    tmp_data = [tmp_MeanBefore, tmp_MeanDuring, tmp_MeanAfter];

    [KWrateP(i), ~, KWrateStat{i}] = kruskalwallis(tmp_data, [], 'off');
    

    subplot(1, n_binSize, i);hold on;grid on;box on;
    b = boxchart(tmp_data, 'Notch','on');
    xticklabels({'Before meal', 'During Meal', 'After meal'});
    ylim([-2,3]);ylabel('Normalized Firing Rate');
    title(['Kurskal-Wallis, P=', num2str(KWrateP(i))]);
end
%%
% lme1 = fitlme(t, ...
%     'MeanFiringDiff_15 ~  1 + (1|RatID) + (1|ExpDay)+ (1|ExpDay:RatID)')
% lme2 = fitlme(t, ...
%     'MeanFiringDiff_15 ~ 1 + RatID*ExpDay + (1|CellID)')
% tt_b = table(t.MeanFiringBefore_15, t.CellID, t.RatID, t.ExpDay, ones(N_Cells, 1), ...
%     'VariableNames',{'MeanFiring', 'CellID', 'RatID', 'ExpDay', 'Eaten'});
% tt_a = table(t.MeanFiringAfter_15, t.CellID, t.RatID, t.ExpDay, 2*ones(N_Cells, 1), ...
%     'VariableNames',{'MeanFiring', 'CellID', 'RatID', 'ExpDay', 'Eaten'});
% tt = [tt_b; tt_a];
% tt.Eaten = categorical(tt.Eaten);
% lme1 = fitlme(tt, ...
%     'MeanFiring ~ Eaten + (1|RatID) + (1|ExpDay)+ (1|ExpDay:RatID)');
% lme2 = fitlme(tt, ...
%     'MeanFiring ~ Eaten*RatID*ExpDay');
%%
KWslopeP = zeros(1, n_binSize);
KWslopeStat = cell(1, n_binSize);
fig4 = figure('Name','FiringSlopeAroundMeal', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    tmp_SlopeBefore = t.(['SlopeFiringBefore_', num2str(binSize_vec(i))]);
    tmp_SlopeDuring = t.(['SlopeFiringDuring_', num2str(binSize_vec(i))]);
    tmp_SlopeAfter = t.(['SlopeFiringAfter_', num2str(binSize_vec(i))]);
    
    tmp_data = [tmp_SlopeBefore, tmp_SlopeDuring, tmp_SlopeAfter]; 

    [KWslopeP(i), ~, KWslopeStat{i}] = kruskalwallis(tmp_data, [], 'off');


    subplot(1, n_binSize, i);hold on;grid on;box on;
    boxchart(tmp_data, 'Notch','on');
    xticklabels({'Before meal', 'During meal', 'After meal'});
    ylim([-6,8]);ylabel('Firing Rate slope');
    title(['Kurskal-Wallis, P=', num2str(KWslopeP(i))]);

end
%%
fig4 = figure('Name','FiringSlopeAroundMeal', 'OuterPosition',[50,50,900,1200]);

for i = 1:n_binSize
    tmp_FitBefore = [t.(['FittedFiringBefore_', num2str(binSize_vec(i))]){:}]';
    tmp_FitBeforeMean = mean(tmp_FitBefore);
    tmp_FitBeforeSem = std(tmp_FitBefore)./sqrt(N_Cells);

    tmp_FitDuring = [t.(['FittedFiringDuring_', num2str(binSize_vec(i))]){:}]';
    tmp_FitDuringMean = mean(tmp_FitDuring);
    tmp_FitDuringSem = std(tmp_FitDuring)./sqrt(N_Cells);

    tmp_FitAfter = [t.(['FittedFiringAfter_', num2str(binSize_vec(i))]){:}]';
    tmp_FitAfterMean = mean(tmp_FitAfter);
    tmp_FitAfterSem = std(tmp_FitAfter)./sqrt(N_Cells);


    ax(i*3-2) = subplot(n_binSize, 3, i*3-2);hold on;grid on;box on;
    shadedErrorBar(tBeforeVec{i}, tmp_FitBeforeMean, tmp_FitBeforeSem);
    ax(i*3-1) = subplot(n_binSize, 3, i*3-1);hold on;grid on;box on;
    shadedErrorBar(tDuringVec{i}, tmp_FitDuringMean, tmp_FitDuringSem);
    ax(i*3) = subplot(n_binSize, 3, i*3);hold on;grid on;box on;
    shadedErrorBar(tAfterVec{i}, tmp_FitAfterMean, tmp_FitAfterSem);

end
linkaxes(ax, 'y')
%%
fige1 = figure('Name','Fig_FiringAroundMealExamples', 'OuterPosition',[50,50,600,1200]);
i = 4;

subplot(4,1,1);
    tmp_z = mean(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1);
    tmp_sem = std(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1)./sqrt(N_Cells);

    shadedErrorBar(tMealVec{i}, tmp_z, tmp_sem);hold on;grid on;box on;
    xregion(-1, 0, 'FaceColor','y');
    xregion(0, 1, 'FaceColor','c');
    xregion(1, 2, 'FaceColor','g');

    legend({'Mean & Standard Error', 'Before Meal', 'During Meal', 'After Meal'}, ...
        'Location','south', 'Orientation','horizontal');

    xlim([-4, 6]);xlabel('Time');
    ylim([-1, 1]);ylabel('Normalized Firing Rate (z-score)');
    title('Average normalized firing of neurons around first meal');

subplot(4,1,2:4);
    tmp_MeanFiringBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    [~, tmp_order] = sort(tmp_MeanFiringBefore, 'ascend');
    tmp_FiringMeal = vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:});
    tmp_FiringMealSorted = tmp_FiringMeal(tmp_order,:);

    imagesc(tMealVec{i}, 1:N_Cells, tmp_FiringMealSorted);hold on;box on;
    % set(gca, 'Ydir','normal');
    xlim([-4, 6]);xlabel('Time');
    ylim([1, N_Cells]);ylabel('Cells');
    title('Ordered normalized firing of singlue neurons around first meal');
    colormap jet;
    c = colorbar('southoutside');
    c.Label.String = 'Normalzied Firing Rate (z-score)';

% printjpg(fige1);printeps(fige1);

%%
fige2 = figure('Name','Fig_FiringAroundMealBox', 'OuterPosition',[50,50,400,800]);
subplot(2,1,1)
    tmp_MeanBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    tmp_MeanDuring = t.(['MeanFiringDuring_', num2str(binSize_vec(i))]);
    tmp_MeanAfter = t.(['MeanFiringAfter_', num2str(binSize_vec(i))]);
    tmp_dataMean = [tmp_MeanBefore, tmp_MeanDuring, tmp_MeanAfter];

    boxchart(tmp_dataMean, 'Notch','on');hold on;grid on;box on;
    ylabel('Normalized Firing Rate (z-score)');ylim([-2,4]);
    xticklabels({'Before meal', 'During meal', 'After meal'});

subplot(2,1,2)
    tmp_SlopeBefore = t.(['SlopeFiringBefore_', num2str(binSize_vec(i))]);
    tmp_SlopeDuring = t.(['SlopeFiringDuring_', num2str(binSize_vec(i))]);
    tmp_SlopeAfter = t.(['SlopeFiringAfter_', num2str(binSize_vec(i))]);
    tmp_dataSlope = [tmp_SlopeBefore, tmp_SlopeDuring, tmp_SlopeAfter];

    boxchart(tmp_dataSlope, 'Notch','on');hold on;grid on;box on;
    ylabel('Slope');ylim([-7,7]);
    xticklabels({'Before meal', 'During meal', 'After meal'});
 
% printjpg(fige2);printeps(fige2);

%%
fige3 = figure('Name','Fig_FiringFit', 'OuterPosition',[50,50,400,400]);
    tmp_FitBefore = [t.(['FittedFiringBefore_', num2str(binSize_vec(i))]){:}]';
    tmp_FitBeforeMean = mean(tmp_FitBefore);
    tmp_FitBeforeSem = std(tmp_FitBefore)./sqrt(N_Cells);

    tmp_FitDuring = [t.(['FittedFiringDuring_', num2str(binSize_vec(i))]){:}]';
    tmp_FitDuringMean = mean(tmp_FitDuring);
    tmp_FitDuringSem = std(tmp_FitDuring)./sqrt(N_Cells);

    tmp_FitAfter = [t.(['FittedFiringAfter_', num2str(binSize_vec(i))]){:}]';
    tmp_FitAfterMean = mean(tmp_FitAfter);
    tmp_FitAfterSem = std(tmp_FitAfter)./sqrt(N_Cells);

    ax(1) = subplot(1,3,1);hold on;grid on;box on;
    shadedErrorBar(tBeforeVec{i}, tmp_FitBeforeMean, tmp_FitBeforeSem);
    ax(2) = subplot(1,3,2);hold on;grid on;box on;
    shadedErrorBar(tDuringVec{i}, tmp_FitDuringMean, tmp_FitDuringSem);
    ax(3) = subplot(1,3,3);hold on;grid on;box on;
    shadedErrorBar(tAfterVec{i}, tmp_FitAfterMean, tmp_FitAfterSem);

    ylim([-0.5,0.7]);

    linkaxes(ax, 'y');

    % printjpg(fige3);printeps(fige3);
%%
save('Gut_ENS_ChronicExp2_FiringRateMealSlope.mat', 'T');

