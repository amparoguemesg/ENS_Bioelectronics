%% Housekeeping
close all;clearvars;clc;

%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
f = [34];
f_length = length(f);
for idx = 1:f_length
    idxFile = f(idx);
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);

    % KiloSort Compatible Channel Map
    Nchannels = recordingInfo.intanChNum;
    chanMap0ind = [recordingInfo.nrscpChMap{:}]';
    chanMap = chanMap0ind + 1;
    connected = true(Nchannels, 1);
    fs = recordingInfo.recordingSamplingRate;
    kcoords = [recordingInfo.nrscpChMapK{:}]';
    xcoords = [recordingInfo.nrscpChMapX{:}]';
    ycoords = [recordingInfo.nrscpChMapY{:}]';

    save([nameFile, '_KS_ChMap.mat'], ...
        'Nchannels', 'chanMap', 'chanMap0ind', 'connected', 'fs', ...
        'kcoords', 'xcoords', 'ycoords');
    
    dummy = 1;
    cd ..
end
