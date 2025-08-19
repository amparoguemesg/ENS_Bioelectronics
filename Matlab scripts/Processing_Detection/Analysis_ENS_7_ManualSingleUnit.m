%% Housekeeping
close all;clearvars;clc;
%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
tbl = readtable('Gut_ENS_Progress.xlsx');
%%
f = [25, 27, 29, 32, 34, 36];
f_length = length(f);
for idx = 1:f_length
    idxFile = f(idx);
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);
    
    % Load KiloSort results
    dir_KS_name = [pwd, '/KS_OutputMed_', nameFile];
    sp = loadKSdir(dir_KS_name);   
    sp_dpt = readNPY([pwd, '/KS_OutputMed_', nameFile, '/spike_times.npy']);
    pc_features = readNPY([pwd, '/KS_OutputMed_', nameFile, '/pc_features.npy']);
    pc_features_ind = readNPY([pwd, '/KS_OutputMed_', nameFile, '/pc_feature_ind.npy']);
    load([pwd, '/', nameFile, '_KS_ChMap.mat']);
    N_chs = Nchannels;

    chanMap0ind_connected = chanMap0ind(connected);
    [N_clusters, N_neighborChans] = size(pc_features_ind);
    pc_feature_chs0ind = zeros(N_clusters, N_neighborChans);
    for idxCluster = 1:N_clusters
        pc_feature_chs0ind(idxCluster,:) = ...
            chanMap0ind_connected(pc_features_ind(idxCluster,:)+1);
    end
    
    % Manually selected clusters
    su_ids = str2num(string(tbl.KiloSort{strcmp(nameFile, tbl.RecordingName)}))';
    N_cells = numel(su_ids);

    sp_clu = sp.clu;
    sp_time = sp.st;
    su_n = zeros(N_cells, 1);
    su_t = cell(N_cells, 1);
    su_dpt = cell(N_cells, 1);
    su_pc = cell(N_cells, 1);
    for idxCell = 1:N_cells
        tmp_idx = sp_clu == su_ids(idxCell);
        su_n(idxCell) = sum(tmp_idx);
        su_t{idxCell} = sp_time(tmp_idx);
        su_dpt{idxCell} = sp_dpt(tmp_idx);
        su_pc{idxCell} = pc_features(tmp_idx,:,:);
    end
    su_pc_chs = pc_feature_chs0ind(su_ids+1, :);

    % Remove fast artifacts
    dp_refract = 30;
    for idxCell = 1:N_cells
        tmp_idx = [0; diff(su_dpt{idxCell})] >= dp_refract;
        su_n(idxCell) = sum(tmp_idx);
        su_t{idxCell} = su_t{idxCell}(tmp_idx);
        su_dpt{idxCell} = su_dpt{idxCell}(tmp_idx);
        su_pc{idxCell} = su_pc{idxCell}(tmp_idx,:,:);
    end
    
    % calculate ACG & CCG
    ACG = struct;
    binSize = 1e-3;
    nBins = 51;
    convWin = 11;
    alpha = 0.05;
    for i = 1:N_cells
        for j = 1:N_cells
            [tmp_c, tmp_b] = fma_CCG( vertcat(su_t{i}, su_t{j}), ...
                vertcat(ones(size(su_t{i})), 2*ones(size(su_t{j}))), ...
                'binSize',binSize, 'nBins',nBins);
            tmp_b = tmp_b * 1e3;
            if i == j
                tmp_c = tmp_c(:,1,1);
            else
                tmp_c = tmp_c(:,1,2);
            end
            [~, tmp_pred, ~] = cch_conv(round(tmp_c), convWin);

            tmp_hi = poissinv(1-alpha, tmp_pred);
            tmp_lo = poissinv(alpha, tmp_pred);

            ACG(i,j).Cell1 = su_ids(i);
            ACG(i,j).Cell2 = su_ids(j);
            ACG(i,j).c = tmp_c;
            ACG(i,j).b = tmp_b;
            ACG(i,j).lo = tmp_lo;
            ACG(i,j).hi = tmp_hi;
        end
    end

    % Save
    save([nameFile, '_ManualSingleUnit.mat'], ...
        'su*', 'ACG', 'binSize', 'nBins', 'convWin', 'alpha');

    dummy = 1;
    cd ..
end
