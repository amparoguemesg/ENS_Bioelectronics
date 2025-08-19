%% Housekeeping
close all;clearvars;clc;
%%
load("Analysis_ENS_Data.mat");
% remove Urethane_distension
dataInfo(23) = [];
N = N - 1;
%%
Thresholds = [4, 7];
%%
f = [34];
f_length = length(f);
for idx = 1:f_length
    idxFile = f(idx);
    nameFile = dataInfo(idxFile).recordingName;
    cd(nameFile);
    load([nameFile, '_recordingInfo.mat']);
    load([nameFile, '_NoiseFloor.mat']);
    load([nameFile, '_artifacts.mat']);

    [PeaksAll, VoltagesAll] = fil2MUA_Dev(recordingInfo.medName, ...
        recordingInfo.recordingSamplingRate, recordingInfo.intanChNum, ...
        recordingInfo.mtlbChList, mean(noiseFloor_median,2), ...
        Thresholds(1), Thresholds(2));

    PeaksQuiet = cell(recordingInfo.intanChNum, 1);
    VoltagesQuiet = cell(size(PeaksQuiet));
    for idxCh = 1:recordingInfo.intanChNum
        [PeaksQuiet{idxCh}, in_state_idx] = restrictMUAState(PeaksAll{idxCh}, ...
            intervalQuiet_dp);
        VoltagesQuiet{idxCh} = VoltagesAll{idxCh}(in_state_idx);
    end
    save([nameFile, '_MUADetect.mat'], ...
        'PeaksQuiet', 'VoltagesQuiet', 'Thresholds');

    clu = 1;
    for idxCh = 1:recordingInfo.intanChNum
        if size(PeaksQuiet{idxCh},1)==1
            resQuiet = PeaksQuiet{idxCh}';
        else
            resQuiet = PeaksQuiet{idxCh};
        end
        cluQuiet = vertcat(clu, ...
            recordingInfo.mtlbChList(idxCh).*ones(length(resQuiet),1));
    
        resName = [nameFile, '_CH', num2str(recordingInfo.mtlbChList(idxCh)), ...
            '_MUA', '.res'];
        dlmwrite(resName, resQuiet, 'precision',10);

        cluName = [nameFile, '_CH', num2str(recordingInfo.mtlbChList(idxCh)), ...
            '_MUA', '.clu'];
        dlmwrite(cluName, cluQuiet, 'precision',10);
    end
    
    nameFolderEvt = ['Events_', nameFile, '_MUA'];
    if isempty(dir(nameFolderEvt)); mkdir(nameFolderEvt); end
    movefile('*_MUA.res', nameFolderEvt);
    movefile('*_MUA.clu', nameFolderEvt);
    
    dummy = 1;
    cd ..
end