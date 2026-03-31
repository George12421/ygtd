% ---------------- generateTestData.m (双变量多场景测试集版) ----------------
% 用于生成验证算法精度的 Ground Truth 测试数据集 (树状级联序列)
% 适配负载与新能源双重波动的鲁棒性测试
% ----------------------------------------------------

clear
params = getParams_v3(true); 
flowType = params.flowType;
mpc = params.caseMatpowerStuct;

% --- 测试集核心参数：定义测试场景矩阵 ---
% 每一行代表一个测试场景: [负载比例 LoadScale, 新能源比例 RenewScale]
test_scenarios = [
    1.00, 0.00;  % 场景 1: 纯基准情况 (无新能源)
    1.10, 0.20;  % 场景 2: 高负载 + 低新能源 (系统高压状态)
    0.90, 0.80;  % 场景 3: 低负载 + 高新能源 (易发潮流倒送)
    1.05, 0.50   % 场景 4: 中等负载 + 中等新能源
];
num_scenarios = size(test_scenarios, 1);
max_depth = 3; % 级联最大深度  14-7，39-5，118-3

% --- 自动适配新能源节点 (与训练数据生成一致) ---
all_gen_buses = mpc.gen(:, 1);
solar_template = struct();
wind_template = struct();

if strcmp(params.caseName, 'IEEE14')
    solar_template.bus_numbers = [2]; 
    wind_template.bus_numbers  = [3];
    P_solar_max_list = [40]; 
    P_wind_max_list  = [60];
elseif strcmp(params.caseName, 'IEEE39')
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

renew_buses = [solar_template.bus_numbers, wind_template.bus_numbers];
solar_template.conventional_gens = find(~ismember(all_gen_buses, renew_buses));

% 准备基准电网
[mpc_base, ~, ~] = prepareGridStruct(mpc, flowType); 
numFailures = params.numSamples;
initialFailureIndicesAll = getInitialFailures(numFailures, mpc_base.numInitialLines, 1);

% 基础保存目录
baseTestDataDir = fullfile(params.dataPath, flowType, params.caseName, 'test_cascade_data');
fprintf('信息: 准备生成 Ground Truth 测试序列 (最大级联深度 M=%d)。\n', max_depth);

% 限制并行worker数以避免内存不足
numWorkers = 8; % 根据您的CPU和内存调整
if isempty(gcp('nocreate'))
    parpool(numWorkers);
end

% =========================================================================
% 外层循环：遍历设定的测试场景
% =========================================================================
for s_idx = 1:num_scenarios
    curr_L = test_scenarios(s_idx, 1);
    curr_R = test_scenarios(s_idx, 2);
    
    fprintf('\n======================================================\n');
    fprintf('>>> 开始生成测试场景 %d/%d: 负载 %.2f, 新能源 %.2f <<<\n', s_idx, num_scenarios, curr_L, curr_R);
    fprintf('======================================================\n');
    
    % A. 构建当前场景的新能源状态
    temp_solar = solar_template;
    temp_solar.power_mw_list = P_solar_max_list * curr_R;
    temp_wind = wind_template;
    temp_wind.power_mw_list = P_wind_max_list * curr_R;
    
    % B. 应用新能源和负载
    mpc_with_re = applyHybridRenewables(mpc_base, temp_solar, temp_wind);
    mpc_final = scaleLoadGenProfile(mpc_with_re, curr_L);
    
    % C. 为当前场景创建独立子文件夹
    l_str = strrep(sprintf('%.2f', curr_L), '.', 'p');
    r_str = strrep(sprintf('%.2f', curr_R), '.', 'p');
    scenarioDirName = sprintf('scenario_L%s_R%s', l_str, r_str);
    testDataDir = fullfile(baseTestDataDir, scenarioDirName);
    
    if ~exist(testDataDir, 'dir'); mkdir(testDataDir); end
    total_sequences_generated = 0;
    
    % =========================================================================
    % 内层循环：遍历 N-1 初始故障
    % =========================================================================
    for jj_fail = 0:numFailures-1
        kk_fail = jj_fail; 
        
        initialFailureIndices = initialFailureIndicesAll(:,jj_fail+1);
        
        % 生成树状级联序列
        all_sequences = powerCascadeTree(mpc_final, initialFailureIndices, flowType, max_depth);
        
        if ~isempty(all_sequences)
            failureDir = fullfile(testDataDir, sprintf('failure_%04d', kk_fail));
            if ~exist(failureDir, 'dir'); mkdir(failureDir); end
            
            for seq_idx = 1:length(all_sequences)
                fileName = sprintf('cascade_seq_%04d.mat', seq_idx);
                saveLocation = fullfile(failureDir, fileName);
                saveData(saveLocation, all_sequences{seq_idx});
            end
            
            total_sequences_generated = total_sequences_generated + length(all_sequences);
            fprintf('  [线路 %04d] 生成了 %d 条级联分支。\n', kk_fail, length(all_sequences));
        end
    end
    fprintf('--- 场景 %d 测试数据生成完毕！共生成 %d 条 Ground Truth 序列 ---\n', s_idx, total_sequences_generated);
end

disp(' ');
disp('--- 所有多场景测试数据集生成完毕！ ---');