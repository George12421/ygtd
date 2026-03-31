% ---------------- generateData_LHS_v4.m (双重波动 + LHS 采样版) ----------------
% 结构: failure_XXXX / data_L[load]_R[renew].mat
% 适用: IEEE 14, 39, 118 系统的级联故障数据生成

clear
% 1. 读取仿真参数
params         = getParams_v3(true); 
flowType       = params.flowType;
mpc            = params.caseMatpowerStuct;
numFailures    = params.numSamples; % 故障数
enableSaveData = params.enableSaveData;

% =========================================================================
% 核心设置 1：LHS 联合采样 (负载 + 新能源)
% =========================================================================
numJointSamples = 5000; % 建议样本量，确保 ICA 训练精度 [cite: 5, 203]
loadBounds  = [0.90, 1.10];  % 负载波动范围
renewBounds = [0.00, 1.00];  % 新能源出力范围

fprintf('信息: 正在生成 LHS 采样矩阵 (%d 样本)...\n', numJointSamples);
% 使用拉丁超立方采样，保证负载与新能源波动的统计独立性 [cite: 78, 200]
X_lhs = lhsdesign(numJointSamples, 2); 
loadSamples  = loadBounds(1) + (loadBounds(2) - loadBounds(1)) * X_lhs(:,1);
renewSamples = renewBounds(1) + (renewBounds(2) - renewBounds(1)) * X_lhs(:,2);

% =========================================================================
% 核心设置 2：自动适配系统节点 (修正 IEEE 39 错误)
% =========================================================================
all_gen_buses = mpc.gen(:, 1);
solar_template = struct();
wind_template = struct();

if strcmp(params.caseName, 'IEEE14')
    solar_template.bus_numbers = [2]; 
    wind_template.bus_numbers  = [3];
    P_solar_max_list = [40]; 
    P_wind_max_list  = [60];
elseif strcmp(params.caseName, 'IEEE39')
    % 修正：IEEE 39 的发电机母线位于 30-39 之间
    solar_template.bus_numbers = [30, 31]; 
    wind_template.bus_numbers  = [32, 33];
    P_solar_max_list = [300, 300];
    P_wind_max_list  = [400, 400];
elseif strcmp(params.caseName, 'IEEE118')
    solar_template.bus_numbers = [10, 12, 25]; 
    wind_template.bus_numbers  = [26, 31, 46];
    P_solar_max_list = [150, 150, 150]; 
    P_wind_max_list  = [200, 200, 200];
else
    error('请配置当前系统 %s 的新能源母线', params.caseName);
end

% 自动计算剩余常规机组索引用于功率平衡
renew_buses = [solar_template.bus_numbers, wind_template.bus_numbers];
solar_template.conventional_gens = find(~ismember(all_gen_buses, renew_buses));

% -----------------------------------------------------------------
% 准备基础网架
[mpc_base, ~, ~] = prepareGridStruct(mpc, flowType); 
initialFailureIndicesAll = getInitialFailures(numFailures, mpc_base.numInitialLines, 1);

% 获取全局保存路径
[sampleDataDir, ~] = getSaveDirectory(params, 1.0);
caseRootDirectory = fileparts(fileparts(sampleDataDir));
globalLinksetDirectory = fullfile(caseRootDirectory, 'linkset');
if ~exist(globalLinksetDirectory, 'dir'); mkdir(globalLinksetDirectory); end
saveLinksetAndCapacity(mpc_base, globalLinksetDirectory, 0, params.caseName);

% =========================================================================
% 主循环：Faults -> Parfor (Load & Renewable LHS Samples)
% =========================================================================
% 限制并行worker数以避免内存不足
numWorkers = 12; % 根据您的CPU和内存调整
if isempty(gcp('nocreate'))
    parpool(numWorkers);
end

for jj_fail = 0:numFailures-1
    kk_fail = jj_fail;
    fprintf('\n--- 开始处理初始故障: %04d (LHS 并行采样开始) ---\n', kk_fail);
    
    initialFailureIndices = initialFailureIndicesAll(:, jj_fail+1);
    failureString = sprintf('failure_%04d', kk_fail);
    dataDirectory = fullfile(caseRootDirectory, failureString);
    if ~exist(dataDirectory, 'dir'); mkdir(dataDirectory); end
    
    totalTime = 0; numFails = 0;

    % 并行火力全开处理联合采样点
    parfor ii = 1:numJointSamples
        curr_L = loadSamples(ii);
        curr_R = renewSamples(ii);
        
        % A. 构建当前样本新能源出力
        temp_solar = solar_template;
        temp_solar.power_mw_list = P_solar_max_list * curr_R;
        temp_wind = wind_template;
        temp_wind.power_mw_list = P_wind_max_list * curr_R;
        
        % B. 先应用新能源，再应用负载缩放
        mpc_with_re = applyHybridRenewables(mpc_base, temp_solar, temp_wind);
        mpc_final = scaleLoadGenProfile(mpc_with_re, curr_L);
        
        % C. 文件名包含双重参数，方便提取
        l_str = strrep(sprintf('%.4f', curr_L), '.', 'p');
        r_str = strrep(sprintf('%.4f', curr_R), '.', 'p');
        fileName = ['data_L' l_str '_R' r_str '.mat'];
        saveLocation = fullfile(dataDirectory, fileName);
        
        % D. 运行仿真
        t_start = tic;
        success = powerCascade(mpc_final, initialFailureIndices, saveLocation, enableSaveData, flowType);
        t_end = toc(t_start);
        
        totalTime = totalTime + t_end;
        if success ~= 1, numFails = numFails + 1; end
    end

    fprintf('   完成故障 %d: 平均耗时 %.2f ms, 失败 %d\n', ...
            kk_fail, (1000 * totalTime / numJointSamples), numFails);
end


