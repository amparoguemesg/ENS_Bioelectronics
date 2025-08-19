%% Housekeeping
close all;clearvars;clc;

%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
load('Analysis_ENS_ChronicExp2MealTimes.mat');
%%
binSize = 30;
binEdges = 0:binSize:30*60;
binX = (binEdges(1:end-1) + binEdges(2:end))./2;
binX = binX./60;
%%
for idxFile = 10:21
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);
    load([nameFile, '_artifacts.mat'], 'intervalArtifacts');
    intervalArtifactsMin = intervalArtifacts ./ 60;

    MealTimeMin = tblMealTime.MealTimeMinutes{ ...
        strcmp(nameFile, tblMealTime.RecordingName)};
    intervalMealTimeMin = [MealTimeMin; MealTimeMin+1]';

    % Load single units
    load([nameFile, '_ManualSingleUnit.mat'], 'su_t');    
    nCells = length(su_t);
    su_th = cell(nCells,1);
    
    for idxCell = 1:nCells
        su_th{idxCell} = histcounts(su_t{idxCell}, binEdges)./binSize;
    end
    
    fig = figure('Name',['Figure_', nameFile, '_FiringHistogram'], ...
        'OuterPosition',[50,50,2400,600]);
    for idxCell = 1:nCells
        plot(binX, su_th{idxCell});hold on;grid on;
    end
    xlim([0, 30]);
    xregion(intervalArtifactsMin, 'FaceColor','r');
    xregion(intervalMealTimeMin, 'FaceColor','b');
    xlabel('Time (min)');ylabel('Instantaneous Firing Rate(Hz)');
    title(nameFile, 'Interpreter','none');

    printjpg(fig);
    close all;
    cd ..
end
