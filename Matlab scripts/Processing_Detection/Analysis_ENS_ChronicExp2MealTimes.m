%% Housekeeping
close all;clearvars;clc;
%% 4 rats on 3 days -- first received from London
MealTimeMins = cell(4,3);

% day 1
MealTimeMins{1,1} = [15, 16, 28];
MealTimeMins{2,1} = [15, 16];
MealTimeMins{3,1} = [15 ,21, 22, 23];
MealTimeMins{4,1} = [20, 21];

% day 8
MealTimeMins{1,2} = [15, 16, 28];
MealTimeMins{2,2} = [15, 16];
MealTimeMins{3,2} = [15 ,21, 22, 23];
MealTimeMins{4,2} = [20, 21];

% day 12
MealTimeMins{1,3} = [15, 16, 21, 22];
MealTimeMins{2,3} = [20, 28];
MealTimeMins{3,3} = [15, 16, 17];
MealTimeMins{4,3} = [15, 16];
%% Updated from london
MealTimeExactSec = cell(4,3);
% day 1
MealTimeExactSec{1,1} = [15*60+15, 28*60+15];
MealTimeExactSec{2,1} = [15*60+45];
MealTimeExactSec{3,1} = [15*60+48];
MealTimeExactSec{4,1} = [20*60+10];
% day 8
MealTimeExactSec{1,2} = [15*60+5];
MealTimeExactSec{2,2} = [15*60+17];
MealTimeExactSec{3,2} = [15*60+5];
MealTimeExactSec{4,2} = [15*60+47];
% day 12
MealTimeExactSec{1,3} = [15*60+9, 21*60+32];
MealTimeExactSec{2,3} = [20*60+9, 28*60+7];
MealTimeExactSec{3,3} = [15*60+23, 15*60+58];
MealTimeExactSec{4,3} = [15*60+21];

%%
RecordingName = cell(12,1);
RecordingName{1} = 'Chronic_experiment2_day1_rat1';
RecordingName{2} = 'Chronic_experiment2_day1_rat2';
RecordingName{3} = 'Chronic_experiment2_day1_rat3';
RecordingName{4} = 'Chronic_experiment2_day1_rat4';
RecordingName{5} = 'Chronic_experiment2_day8_rat1';
RecordingName{6} = 'Chronic_experiment2_day8_rat2';
RecordingName{7} = 'Chronic_experiment2_day8_rat3';
RecordingName{8} = 'Chronic_experiment2_day8_rat4';
RecordingName{9} = 'Chronic_experiment2_day12_rat1';
RecordingName{10} = 'Chronic_experiment2_day12_rat2';
RecordingName{11} = 'Chronic_experiment2_day12_rat3';
RecordingName{12} = 'Chronic_experiment2_day12_rat4';

MealTimeMinutes = MealTimeMins(:);

MealTimeExactSec = MealTimeExactSec(:);

tblMealTime = table(RecordingName, MealTimeMinutes, MealTimeExactSec);

save('Analysis_ENS_ChronicExp2MealTimes.mat', 'tblMealTime');