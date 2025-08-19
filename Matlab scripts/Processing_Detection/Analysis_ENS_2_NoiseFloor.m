%% Housekeeping
close all;clearvars;clc;
%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%% Define segments of recording for noise floor calculation
Rs = 1e4;
%%
f = [32, 34, 36];
f_length = length(f);
for idx = 1:f_length
    idxFile = (f(idx));
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);

    % manually scored artifact intervals
    dir_art_evt = dir('*artifacts.art.evt');
    art_struct = LoadEvents(dir_art_evt.name);
    numArtInts = length(art_struct.time) / 2;
    intervalArtifacts = reshape(art_struct.time, 2, numArtInts);
    intervalArtifacts = intervalArtifacts';
    intervalArtifacts_dp = round(intervalArtifacts.* Rs);
    labelArtiacts = evt2label(intervalArtifacts_dp, ...
        recordingInfo.datLength);
    labelQuiet = ~labelArtiacts;
    intervalQuiet_dp = label2evt(labelQuiet);
    intervalQuiet = intervalQuiet_dp./ Rs;
    save([nameFile, '_artifacts.mat'], ...
        'labelArtiacts', 'intervalArtifacts_dp', 'intervalArtifacts', ...
        'labelQuiet', 'intervalQuiet_dp', 'intervalQuiet');

    % estimate noise floor from long state intervals
    intervalQuietDuration = intervalQuiet_dp(:,2) - intervalQuiet_dp(:,1);
    [~, interval_idx] = sort(intervalQuietDuration, 'descend');
    N_nf_interval = 5;
    nf_interval = intervalQuiet_dp(interval_idx(1:N_nf_interval),:);
    nf_interval_dur = intervalQuietDuration(interval_idx(1:N_nf_interval));
    nf_interval_mid = round((nf_interval(:,1) + nf_interval(:,2))/2);

    noiseFloor_median = zeros(recordingInfo.intanChNum, N_nf_interval);
    noiseFloor_mean = zeros(size(noiseFloor_median));
    for j = 1:N_nf_interval
        tmp_data = Dat_tracker(recordingInfo.medName, nf_interval_mid(j), ...
            round(nf_interval_dur(j) * 0.95), recordingInfo.intanChNum);
        noiseFloor_median(:,j) = median(abs(tmp_data),2)./0.6745;
        noiseFloor_mean(:,j) = mean(abs(tmp_data),2)./0.6745;
    end

    save([nameFile, '_noiseFloor.mat'], ...
        'noiseFloor*', 'nf_interval*');

    
    dummy = 1;
    cd ..
end