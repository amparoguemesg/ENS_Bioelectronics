%% Housekeeping
close all;clearvars;clc;
%%
dir_info = dir("*_recordingInfo.mat");
N_recordings = length(dir_info);
load('../Analysis_ENS_ChronicExp2MealTimes.mat');
%%
n_cells = zeros(1, N_recordings);

binSize_vec = [10, 20, 30, 40, 50, 60];
n_binSize = length(binSize_vec);

RatID = [];
ExpDay = [];
histAll = [];
histAllBefore = [];
histAllAfter = [];
mdlAllBefore = [];
mdlAllAfter = [];
slpAllBefore = [];
slpAllAfter = [];
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
    histFileBefore = cell(nCell, n_binSize);
    histFileAfter = cell(nCell, n_binSize);
    mdlFileBefore = cell(nCell, n_binSize);
    mdlFileAfter = cell(nCell, n_binSize);
    slpFileBefore = zeros(nCell, n_binSize);
    slpFileAfter = zeros(nCell, n_binSize);
    for idxSize = 1:n_binSize
        binEdges = 0:binSize_vec(idxSize):30*60;
        binCenters = (binEdges(1:end-1) + binEdges(2:end))./2;
        binsBefore = binCenters <= FirstMealStart * 60 & ...
            binCenters >= (FirstMealStop - 6) * 60;
        binsAfter = binCenters >= FirstMealStop * 60 & ...
            binCenters <= (FirstMealStop + 5) * 60;
        tFitBefore = linspace(-5, 0, sum(binsBefore));
        tFitAfter = linspace(0, 5, sum(binsAfter));

        for idxCell = 1:nCell
            histFile{idxCell, idxSize} = histcounts(su_t{idxCell}, ...
                binEdges) ./ binSize_vec(idxSize);

            histFileBefore{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsBefore);
            mdlFileBefore{idxCell, idxSize} = fitlm(tFitBefore, ...
                histFileBefore{idxCell, idxSize}, 'linear');
            slpFileBefore(idxCell, idxSize) = ...
                mdlFileBefore{idxCell, idxSize}.Coefficients.Estimate(2);

            histFileAfter{idxCell, idxSize} = histFile{idxCell, idxSize}( ...
                binsAfter);
            mdlFileAfter{idxCell, idxSize} = fitlm(tFitAfter, ...
                histFileAfter{idxCell, idxSize}, 'linear');
            slpFileAfter(idxCell, idxSize) = ...
                mdlFileAfter{idxCell, idxSize}.Coefficients.Estimate(2);

            dummy = 1;
        end
    end
    
    RatID = [RatID; ratIDFile*ones(nCell,1)];
    ExpDay = [ExpDay; expDayFile*ones(nCell,1)];
    histAll = [histAll; histFile];
    histAllBefore = [histAllBefore; histFileBefore];
    histAllAfter = [histAllAfter; histFileAfter];
    mdlAllBefore = [mdlAllBefore; mdlFileBefore];
    mdlAllAfter = [mdlAllAfter; mdlFileAfter];
    slpAllBefore = [slpAllBefore; slpFileBefore];
    slpAllAfter = [slpAllAfter; slpFileAfter];

    dummy = 1;
end
N_Cells = sum(n_cells);
CellID = (1:N_Cells)';

T = table(CellID, RatID, ExpDay);
%%
for idxSize = 1:n_binSize
    T.(['SlopeBeforeMeal_', num2str(binSize_vec(idxSize))]) = ...
        slpAllBefore(:, idxSize);
    T.(['SlopeAfterMeal_', num2str(binSize_vec(idxSize))]) = ...
        slpAllAfter(:, idxSize);
    T.(['LinearFitBeforeMeal_', num2str(binSize_vec(idxSize))]) = ...
        mdlAllBefore(:, idxSize);
    T.(['LinearFitAfterMeal_', num2str(binSize_vec(idxSize))]) = ...
        mdlAllAfter(:, idxSize);
    T.(['FiringBeforeMeal_', num2str(binSize_vec(idxSize))]) = ...
        histAllBefore(:, idxSize);
    T.(['FiringAfterMeal_', num2str(binSize_vec(idxSize))]) = ...
        histAllAfter(:, idxSize);
    T.(['FiringHistogram_', num2str(binSize_vec(idxSize))]) = ...
        histAll(:, idxSize);
end
%%
save('Gut_ENS_ChronicExp2_FiringRateMealSlope.mat', 'T');

