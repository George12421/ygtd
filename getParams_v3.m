function params = getParams_v3(isDataGen)
% define simulation parameters here
dataPath           = './data';
utilsPath          = './utils';
matpowerFolder     = 'D:\bianchengfangzhen\matlab\Matlab 2022b 64bit\toolbox\matpower7.1\matpower7.1';
flowType           = 'dc';
caseName           = 'IEEE118';
numSamples         = 179;             %线路条数            
numInitialFailures = 1;
% loadScaleList      = ["random"];
loadScaleList      = [0.90:0.1:1.10]; 
% loadScaleList      = 1.00;
enableSaveData     = 1;

% pack all parameters into a struct
addpath(utilsPath);
params = defineParamStruct(isDataGen, dataPath, ...
    flowType, caseName, numSamples, numInitialFailures, ...
    loadScaleList, enableSaveData, matpowerFolder);
end