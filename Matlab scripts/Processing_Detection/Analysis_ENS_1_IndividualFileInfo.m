%% Housekeeping
close all;clearvars;clc;
%% ProbeMap
MapIntan = cell(9,1);
MapIntan{1} = [31, 16];
MapIntan{2} = [14, 17, 0, 30];
MapIntan{3} = [13, 18, 1, 29];
MapIntan{4} = [12, 19, 2, 28];
MapIntan{5} = [11, 20, 3, 27];
MapIntan{6} = [10, 21, 4, 26];
MapIntan{7} = [9, 22, 5, 25];
MapIntan{8} = [8, 23, 6, 24];
MapIntan{9} = [7, 15];

MapIntanX = cell(9,1);
MapIntanX{1} = [0, 0];
MapIntanX{2} = [0, 25, 50, 75];
MapIntanX{3} = [0, 25, 50, 75];
MapIntanX{4} = [0, 25, 50, 75];
MapIntanX{5} = [0, 25, 50, 75];
MapIntanX{6} = [0, 25, 50, 75];
MapIntanX{7} = [0, 25, 50, 75];
MapIntanX{8} = [0, 25, 50, 75];
MapIntanX{9} = [0, 0];

MapIntanY = cell(9,1);
MapIntanY{1} = [-20 * 100, -10 * 100];
MapIntanY{2} = [0, 25, 0, 25] + 0 * 200;
MapIntanY{3} = [0, 25, 0, 25] + 1 * 200;
MapIntanY{4} = [0, 25, 0, 25] + 2 * 200;
MapIntanY{5} = [0, 25, 0, 25] + 3 * 200;
MapIntanY{6} = [0, 25, 0, 25] + 4 * 200;
MapIntanY{7} = [0, 25, 0, 25] + 5 * 200;
MapIntanY{8} = [0, 25, 0, 25] + 6 * 200;
MapIntanY{9} = 6 * 200 + [10 * 100, 20 * 100];

MapIntanK = cell(9,1);
MapIntanK{1} = [1, 2];
MapIntanK{2} = [3, 3, 3, 3];
MapIntanK{3} = [4, 4, 4, 4];
MapIntanK{4} = [5, 5, 5, 5];
MapIntanK{5} = [6, 6, 6, 6];
MapIntanK{6} = [7, 7, 7, 7];
MapIntanK{7} = [8, 8, 8, 8];
MapIntanK{8} = [9, 9, 9, 9];
MapIntanK{9} = [10, 11];
%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
f = [25, 27, 29, 32, 34, 36];
f_length = length(f);
for idx = 1:f_length
    idxFile = f(idx);
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    
    recordingInfo = dataInfo(idxFile);
    recordingInfo.mtlbChList = 1:recordingInfo.intanChNum;
    recordingInfo.nrscpChList = recordingInfo.mtlbChList - 1;
    % file specific channel map
    nrscpChMap = cell(9,1);
    nrscpChMapX = cell(9,1);
    nrscpChMapY = cell(9,1);
    nrscpChMapK = cell(9,1);
    for i_row = 1:9
        n_col = length(MapIntan{i_row});
        for i_col = 1:n_col
           nrscpCh = find(recordingInfo.intanChList==MapIntan{i_row}(i_col));
           if ~isempty(nrscpCh)
               nrscpCh = nrscpCh - 1;
               nrscpChMap{i_row} = [nrscpChMap{i_row}, nrscpCh];
               nrscpChMapX{i_row} = [nrscpChMapX{i_row}, ...
                   MapIntanX{i_row}(i_col)];
               nrscpChMapY{i_row} = [nrscpChMapY{i_row}, ...
                   MapIntanY{i_row}(i_col)];
               nrscpChMapK{i_row} = [nrscpChMapK{i_row}, ...
                   MapIntanK{i_row}(i_col)];
           end
        end
    end
    recordingInfo.nrscpChMap = nrscpChMap;
    recordingInfo.nrscpChMapX = nrscpChMapX;
    recordingInfo.nrscpChMapY = nrscpChMapY;
    recordingInfo.nrscpChMapK = nrscpChMapK;
    recordingInfo.MapIntan = MapIntan;
    recordingInfo.MapIntanX = MapIntanX;
    recordingInfo.MapIntanY = MapIntanY;
    recordingInfo.MapIntanK = MapIntanK;
    recordingInfo.Rs_lfp = 1250;
    recordingInfo.datName = dir('*dat.dat').name;
    recordingInfo.lfpName = dir('*lfp.lfp').name;
    recordingInfo.filName = dir('*fil.fil').name;
    recordingInfo.medName = dir('*med.med').name;

    x_med = readmulti(recordingInfo.medName, recordingInfo.intanChNum, 1);
    x_lfp = readmulti(recordingInfo.lfpName, recordingInfo.intanChNum, 1);

    recordingInfo.datLength = length(x_med);
    recordingInfo.lfpLength = length(x_lfp);

    save([nameFile, '_recordingInfo.mat'], "recordingInfo");
    
    
    dummy = 1;

    cd ..
end