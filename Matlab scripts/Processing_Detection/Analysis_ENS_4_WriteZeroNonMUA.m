%% Housekeeping
close all;clearvars;clc;
%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
Rs = 1e4;
spikeWindow = [-8, 8];
%%
f = [34];
f_length = length(f);
for idx = 1:f_length
    idxFile = f(idx);
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);
    load([nameFile, '_PeakFinderRes.mat'], 'chList');
    load([nameFile, '_MUADetect.mat'], 'PeaksQuiet');

    nameMedMUA = [nameFile, '_MedMUA', '.fil'];
    file_length = recordingInfo.datLength;
    file_length_sec = file_length/Rs;

    % read and write
    Peaks = PeaksQuiet;
    n_mua_chs = recordingInfo.intanChNum;
    duration = 20 * Rs;
    s_time = duration / 2;
    for stream = s_time:duration:file_length
        idx_start = stream - s_time + 1;
        idx_stop = idx_start + duration - 1;
        idx = idx_start:idx_stop;

        pos_zero = false(n_mua_chs, length(idx));
        for idx_ch = 1:n_mua_chs
            spks_all = Peaks{idx_ch}';
            spks = spks_all(spks_all >= idx_start & spks_all <= idx_stop);
            if ~isempty(spks)
                spks_win = spks + spikeWindow;
                pos_zero(idx_ch,:) = any( and( idx >= spks_win(:,1), ...
                    idx <= spks_win(:,2) ), 1);
            end
        end
        pos_zero = ~pos_zero;
        data = Dat_tracker(recordingInfo.medName, stream, ...
            duration, n_mua_chs);
        for idx_ch = 1:n_mua_chs
            ch = chList(idx_ch);
            data(ch, pos_zero(idx_ch,:)) = 0;
        end
        mat2dat(data, nameMedMUA);
        disp(cat(2,'there is ', num2str(file_length_sec), ...
            ' seconds of data per channel'));
        disp(cat(2, num2str(floor(1000*(stream/Rs)/file_length_sec)/10), ...
            ' % is done!'));

    end
    residual_length = file_length - (stream + round(duration/2));
    residual_center = round(residual_length/2) + stream + s_time;
    if residual_length > 0
        idx_start = residual_center - round(residual_length/2) + 1;
        idx_stop = idx_start + residual_length - 1;
        idx = idx_start:idx_stop;

        pos_zero = false(n_mua_chs, length(idx));
        for idx_c = 1:n_mua_chs
            spks_all = Peaks{idx_ch}';
            spks = spks_all(spks_all >= idx_start & spks_all <= idx_stop);
            if ~isempty(spks)
                spks_win = spks + spikeWindow;
                pos_zero(idx_ch,:) = any( and( idx >= spks_win(:,1), ...
                    idx <= spks_win(:,2) ), 1);
            end
        end
        pos_zero = ~pos_zero;

        data = Dat_tracker(recordingInfo.medName, residual_center, ...
            residual_length, n_mua_chs);
        for idx_ch = 1:n_mua_chs
            ch = chList(idx_ch);
            data(ch, pos_zero(idx_ch,:)) = 0;
        end
        mat2dat(data', nameMedMUA);
        fprintf('\nConversion complete!\n');
    end


    dummy = 1;
    cd ..
end