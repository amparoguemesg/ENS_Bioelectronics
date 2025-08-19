%% Housekeeping
close all;clearvars;clc;
%%
dir_info = dir("*_recordingInfo.mat");
N_recordings = length(dir_info);
load('../Analysis_ENS_ChronicExp2MealTimes.mat');
%%
n_cells = zeros(1, N_recordings);

binSize_vec = [1, 2, 5, 10, 20];
n_binSize = length(binSize_vec);

RatID = [];
ExpDay = [];

histAll = [];
histMealAll = [];
histBeforeAll = [];
histAfterAll = [];
zHistAll = [];
zHistMealAll = [];
zHistBeforeAll = [];
zHistAfterAll = [];
mzBeforeAll = [];
mzAfterAll = [];
lmBeforeAll = [];
lmAfterAll = [];
slpBeforeAll = [];
slpAfterAll = [];
fitBeforeAll = [];
fitAfterAll = [];

for idxFile = 1:N_recordings

    load(dir_info(idxFile).name);
    nameFile = recordingInfo.recordingName;
    ratIDFile = str2double(nameFile(end));
    expDayFile = str2double(nameFile(24:end-5));

    MealTimeMin = tblMealTime.MealTimeMinutes{ strcmp(nameFile, ...
        tblMealTime.RecordingName)};
    FirstMealStart = MealTimeMin(1);
    FirstMealStop = FirstMealStart + 1;

    load([nameFile, '_ManualSingleUnit.mat'], 'su_t');
    nCell = length(su_t);
    n_cells(idxFile) = nCell;
    
    histFile = cell(nCell, n_binSize);
    histMealFile = cell(nCell, n_binSize);
    histBeforeFile = cell(nCell, n_binSize);
    histAfterFile = cell(nCell, n_binSize);
    zHistFile = cell(nCell, n_binSize);
    zHistMealFile = cell(nCell, n_binSize);
    zHistBeforeFile = cell(nCell, n_binSize);
    zHistAfterFile = cell(nCell, n_binSize);
    mzBeforeFile = zeros(nCell, n_binSize);
    mzAfterFile = zeros(nCell, n_binSize);
    lmBeforeFile = cell(nCell, n_binSize);
    lmAfterFile = cell(nCell, n_binSize);
    slpBeforeFile = zeros(nCell, n_binSize);
    slpAfterFile = zeros(nCell, n_binSize);
    fitBeforeFile = cell(nCell, n_binSize);
    fitAfterFile = cell(nCell, n_binSize);

    tMealVec = cell(n_binSize, 1);
    tBeforeVec = cell(n_binSize, 1);
    tAfterVec = cell(n_binSize, 1);
    for idxSize = 1:n_binSize
        binEdges = 0:binSize_vec(idxSize):30*60;
        binCenters = (binEdges(1:end-1) + binEdges(2:end))./2;

        binsMeal = binCenters >= (FirstMealStop - 8) * 60 & ...
            binCenters <= (FirstMealStop + 8) * 60;
        tMeal = linspace(-8, 8, sum(binsMeal));
        tMealVec{idxSize} = tMeal;
        
        binsBefore = binCenters <= FirstMealStop * 60 & ...
            binCenters >= (FirstMealStop - 1) * 60;
        tBefore = linspace(-1, 0, sum(binsBefore));
        tBeforeVec{idxSize} = tBefore;

        binsAfter = binCenters >= FirstMealStop * 60 & ...
            binCenters <= (FirstMealStop + 1) * 60;
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
    histAfterAll = [histAfterAll; histAfterFile];
    zHistAll = [zHistAll; zHistFile];
    zHistMealAll = [zHistMealAll; zHistMealFile];
    zHistBeforeAll = [zHistBeforeAll; zHistBeforeFile];
    zHistAfterAll = [zHistAfterAll; zHistAfterFile];
    mzBeforeAll = [mzBeforeAll; mzBeforeFile];
    mzAfterAll = [mzAfterAll; mzAfterFile];
    lmBeforeAll = [lmBeforeAll; lmBeforeFile];
    lmAfterAll = [lmAfterAll; lmAfterFile];
    slpBeforeAll = [slpBeforeAll; slpBeforeFile];
    slpAfterAll = [slpAfterAll; slpAfterFile];
    fitBeforeAll = [fitBeforeAll; fitBeforeFile];
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
    T.(['MeanFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        mzAfterAll(:, idxSize);
    T.(['SlopeFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        slpBeforeAll(:, idxSize);
    T.(['SlopeFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        slpAfterAll(:, idxSize);
    T.(['FittedFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        fitBeforeAll(:, idxSize);
    T.(['FittedFiringAfter_', num2str(binSize_vec(idxSize))]) = ...
        fitAfterAll(:, idxSize);
    T.(['LinearModelBefore_', num2str(binSize_vec(idxSize))]) = ...
        lmBeforeAll(:, idxSize);
    T.(['LinearModelAfter_', num2str(binSize_vec(idxSize))]) = ...
        lmAfterAll(:, idxSize);
    T.(['FiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        histBeforeAll(:, idxSize);
    T.(['zFiringBefore_', num2str(binSize_vec(idxSize))]) = ...
        zHistBeforeAll(:, idxSize);
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
N_Cells = height(t);
%%
fig1 = figure('Name','FiringAourndMeal', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    
    tmp_MeanFiringBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    [~, tmp_order] = sort(tmp_MeanFiringBefore, 'descend');
    tmp_FiringMeal = vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:});
    tmp_FiringMealSorted = tmp_FiringMeal(tmp_order,:);

    subplot(1, n_binSize, i);
    imagesc(tMealVec{i}, 1:N_Cells, tmp_FiringMealSorted);hold on;grid on;box on;
    set(gca, 'Ydir','normal');
    xlim([-5, 5]);xlabel('Time');
    ylim([1, N_Cells]);ylabel('Cells');
    title(['Mean Firing BinSize_',num2str(binSize_vec(i)), 's'], ...
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
    xregion(-1, 0, 'FaceColor','g');
    xlim([-5, 5]);xlabel('Time');
    title(['Ordered Firing BinSize_',num2str(binSize_vec(i)), 's'], ...
        'Interpreter','none');
end
%%
statRateCompare = zeros(4, n_binSize);

fig3 = figure('Name','FiringAroundMealBox', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    tmp_MeanBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    tmp_MeanAfter = t.(['MeanFiringAfter_', num2str(binSize_vec(i))]);
    
    [statRateCompare(2,i), statRateCompare(1,i)] = ttest(tmp_MeanBefore, tmp_MeanAfter);
    [statRateCompare(3,i), statRateCompare(4,i)] = signrank(tmp_MeanAfter,tmp_MeanBefore);

    subplot(1, n_binSize, i);hold on;grid on;box on;
    boxchart([tmp_MeanBefore, tmp_MeanAfter], 'Notch','on');
    % scatter(1*ones(N_Cells,1), tmp_MeanBefore, ...
    %     'filled', 'XJitter','randn');
    % scatter(2*ones(N_Cells,1), tmp_MeanAfter, ...
    %     'filled', 'XJitter','randn');
    xticklabels({'During meal', 'After meal'});
end

%%
statSlopeCompare = zeros(4, n_binSize);
statSlopeBefore = zeros(4, n_binSize);
statSlopeAfter = zeros(4, n_binSize);

fig4 = figure('Name','FiringSlopeAroundMeal', 'OuterPosition',[50,50,2400,600]);
for i = 1:n_binSize
    tmp_SlopeBefore = t.(['SlopeFiringBefore_', num2str(binSize_vec(i))]);
    tmp_SlopeAfter = t.(['SlopeFiringAfter_', num2str(binSize_vec(i))]);

    [statSlopeCompare(2,i), statSlopeCompare(1,i)] = ttest(tmp_SlopeBefore, tmp_SlopeAfter);
    [statSlopeCompare(3,i), statSlopeCompare(4,i)] = signrank(tmp_SlopeBefore, tmp_SlopeAfter);

    [statSlopeBefore(2,i), statSlopeBefore(1,i)] = ttest(tmp_SlopeBefore);
    [statSlopeBefore(3,i), statSlopeBefore(4,i)] = signrank(tmp_SlopeBefore);

    subplot(1, n_binSize, i);hold on;grid on;box on;
    boxchart([tmp_SlopeBefore, tmp_SlopeAfter], 'Notch','on');
    xticklabels({'During meal', 'After meal'});

end
%%
fig4 = figure('Name','FiringSlopeAroundMeal', 'OuterPosition',[50,50,600,1200]);

for i = 1:n_binSize
    tmp_FitBefore = [t.(['FittedFiringBefore_', num2str(binSize_vec(i))]){:}]';
    tmp_FitBeforeMean = mean(tmp_FitBefore);
    tmp_FitBeforeSem = std(tmp_FitBefore)./sqrt(N_Cells);
    tmp_FitAfter = [t.(['FittedFiringAfter_', num2str(binSize_vec(i))]){:}]';
    tmp_FitAfterMean = mean(tmp_FitAfter);
    tmp_FitAfterSem = std(tmp_FitAfter)./sqrt(N_Cells);
   
    
    ax(i*2-1) = subplot(n_binSize, 2, i*2-1);hold on;grid on;box on;
    shadedErrorBar(tBeforeVec{i}, tmp_FitBeforeMean, tmp_FitBeforeSem);
    ax(i*2) = subplot(n_binSize, 2, i*2);hold on;grid on;box on;
    shadedErrorBar(tAfterVec{i}, tmp_FitAfterMean, tmp_FitAfterSem);

end
linkaxes(ax, 'y')
%%
fige1 = figure('Name','Fig_FiringAroundMealExamples', 'OuterPosition',[50,50,650,1300]);
i = 5;

subplot(4,1,1);
    tmp_z = mean(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1);
    tmp_sem = std(vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:}),1)./sqrt(N_Cells);

    shadedErrorBar(tMealVec{i}, tmp_z, tmp_sem);hold on;grid on;box on;
    xregion(-1, 0, 'FaceColor','g');

    legend({'Mean & Standard Error', 'Meal Time'}, ...
        'Location','northwest');

    xlim([-5, 5]);xlabel('Time');
    ylim([-1, 1]);ylabel('Normalized Firing Rate (z-score)');
    title('Average normalized firing of neurons around first meal');

subplot(4,1,2:4);
    tmp_MeanFiringBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    [~, tmp_order] = sort(tmp_MeanFiringBefore, 'descend');
    tmp_FiringMeal = vertcat(t.(['zFiringMeal_',num2str(binSize_vec(i))]){:});
    tmp_FiringMealSorted = tmp_FiringMeal(tmp_order,:);

    imagesc(tMealVec{i}, 1:N_Cells, tmp_FiringMealSorted);hold on;box on;
    set(gca, 'Ydir','normal');
    xlim([-5, 5]);xlabel('Time');
    ylim([1, N_Cells]);ylabel('Cells');
    title('Ordered normalized firing of singlue neurons around first meal');
    colormap jet;
    c = colorbar('southoutside');
    c.Label.String = 'Normalzied Firing Rate (z-score)';
    caxis([-1,5]);

fige2 = figure('Name','Fig_FiringAroundMealBox', 'OuterPosition',[50,50,500,500]);
    tmp_MeanBefore = t.(['MeanFiringBefore_', num2str(binSize_vec(i))]);
    tmp_MeanAfter = t.(['MeanFiringAfter_', num2str(binSize_vec(i))]);
    p = signrank(tmp_MeanAfter,tmp_MeanBefore);
    boxchart([tmp_MeanBefore, tmp_MeanAfter], 'Notch','on');hold on;grid on;box on;
    ylabel('Normalized Firing Rate (z-score)');
    xticklabels({'During meal', '1 min after meal'});
    title(['Wilcoxon signed rank test, P=', num2str(p)]);

    printjpg(fige1);printeps(fige1);
    printjpg(fige2);printeps(fige2);
%%
save('Gut_ENS_ChronicExp2_FiringRateMealSlope.mat', 'T');

