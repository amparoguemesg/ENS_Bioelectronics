%% Housekeeping
close all;clearvars;clc;
%%
load('Gut_ENS_ChronicExp2_FiringRateMealSlope.mat');
%% Remove day8_rat3 and day8_rat4
idxRemoval = (T.ExpDay == 8) & (T.RatID == 3 | T.RatID == 4);

t = T(~idxRemoval, :);
t.RatID = categorical(t.RatID);
%%
N = height(t);

binSize_vec = [10, 20, 30, 40, 50, 60];
n_binSize = length(binSize_vec);
%% statistical testing
statBefore = zeros(2, n_binSize);
statAfter = zeros(2, n_binSize);
statCompare = zeros(4, n_binSize);
for i = 1:n_binSize
    [statBefore(2,i), statBefore(1,i)] = ttest( ...
        (t.(['SlopeBeforeMeal_', num2str(binSize_vec(i))])));
    [statAfter(2,i), statAfter(1,i)] = ttest( ...
        (t.(['SlopeAfterMeal_', num2str(binSize_vec(i))])));
    [statCompare(2,i), statCompare(1,i)] = ttest( ...
        t.(['SlopeBeforeMeal_', num2str(binSize_vec(i))]), ...
        t.(['SlopeAfterMeal_', num2str(binSize_vec(i))]));
    [statCompare(3,i), statCompare(4,i)] = signrank( ...
        t.(['SlopeBeforeMeal_', num2str(binSize_vec(i))]), ...
        t.(['SlopeAfterMeal_', num2str(binSize_vec(i))]));

    t.(['SlopeDifference_', num2str(binSize_vec(i))]) = ...
        t.(['SlopeAfterMeal_', num2str(binSize_vec(i))]) - ...
        t.(['SlopeBeforeMeal_', num2str(binSize_vec(i))]);
end
%%
lme = fitlme(t, ...
    'SlopeDifference_10 ~ 1 + (1|RatID) + (1|ExpDay)+ (1|ExpDay:RatID)');
lme = fitlme(t, ...
    'SlopeDifference_10 ~ 1 + (SlopeDifference_10|RatID) + (SlopeDifference_10|ExpDay)+ (SlopeDifference_10|ExpDay:RatID)')
%%

fig1 = figure('Name','FiringRateSlopes', 'OuterPosition',[50,50,2000,400]);
for i = 1:n_binSize
    subplot(1,n_binSize,i);hold on;grid on;box on;
    scatter(1*ones(N,1), t.(['SlopeBeforeMeal_', num2str(binSize_vec(i))]), ...
        'filled', 'XJitter','randn');
    scatter(2*ones(N,1), t.(['SlopeAfterMeal_', num2str(binSize_vec(i))]), ...
        'filled', 'XJitter','randn');
    xticks(1:2);xticklabels({'Before', 'After'});xlim([0,3]);
    ylabel('Slope');ylim([-5,5]);
    title(['BinSize_',num2str(binSize_vec(i)), 's'], 'Interpreter','none');    
end
printjpg(fig1);







