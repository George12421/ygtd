% Filename: applyRenewableScenario.m
% Location: matlabmodels/utils/
% (这是“多节点注入”的升级版)

function mpc_modified = applyRenewableScenario(mpc_original, renewable_scenario)
%% applyRenewableScenario: 将【多个】新能源出力场景应用到MATPOWER case中
%
%   Inputs:
%       mpc_original: 原始的MATPOWER case结构体
%       renewable_scenario: 一个结构体，定义了新能源场景，包含：
%           - .bus_numbers (vector): 新能源接入的【母线编号列表】
%           - .power_mw_list (vector): 对应的【功率列表 (MW)】
%           - .conventional_gens (vector): 作为调节电源的常规发电机索引(在gen矩阵中的行号)
%
%   Outputs:
%       mpc_modified: 应用新能源场景后，修改过的mpc结构体

% 深度复制原始mpc对象，避免修改原始数据
mpc_modified = mpc_original;

% 检查输入的列表长度是否一致
if length(renewable_scenario.bus_numbers) ~= length(renewable_scenario.power_mw_list)
    error('错误: 母线列表和功率列表的长度必须一致！');
end

% --- 步骤 1: 循环注入所有新能源 (作为负负荷) ---
P_new_renewable_total = 0;
for i = 1:length(renewable_scenario.bus_numbers)
    current_bus = renewable_scenario.bus_numbers(i);
    current_power = renewable_scenario.power_mw_list(i);
    
    bus_idx = find(mpc_modified.bus(:, 1) == current_bus);
    
    if isempty(bus_idx)
        error('错误: 在案例中未找到指定的母线编号 %d', current_bus);
    end
    
    % 从母线负荷中减去新能源出力
    original_load_pd = mpc_modified.bus(bus_idx, 3); % PD 在第3列
    mpc_modified.bus(bus_idx, 3) = original_load_pd - current_power;
    
    fprintf('信息: 在母线 %d 注入 %.2f MW 新能源。\n', current_bus, current_power);
    
    % 累加总注入功率
    P_new_renewable_total = P_new_renewable_total + current_power;
end

fprintf('信息: 总共注入 %.2f MW 新能源。\n', P_new_renewable_total);


% --- 步骤 2: 调整常规发电机出力以维持系统总功率平衡 ---
conventional_gens_indices = renewable_scenario.conventional_gens;

% ... (这部分平衡逻辑与之前完全相同) ...

% 检查常规机组是否存在
if any(conventional_gens_indices > size(mpc_modified.gen, 1))
    error('错误: 常规发电机索引超出范围。');
end

% 计算这些常规机组的总出力和总最小出力
total_conventional_power = sum(mpc_modified.gen(conventional_gens_indices, 2)); % PG 在第2列
total_conventional_pmin = sum(mpc_modified.gen(conventional_gens_indices, 10)); % PMIN 在第10列

% 检查是否有足够的下调空间
if (total_conventional_power - P_new_renewable_total) < total_conventional_pmin
    fprintf('警告: 新能源总注入量过大，常规机组下调后将低于其总最小出力。出力将被限制在最小出力。\n');
    P_to_reduce = total_conventional_power - total_conventional_pmin;
else
    P_to_reduce = P_new_renewable_total;
end

% 按比例减少常规机组的出力
power_to_shed_per_unit = P_to_reduce / length(conventional_gens_indices);

reduced_count = 0;
for i = 1:length(conventional_gens_indices)
    idx = conventional_gens_indices(i);
    original_Pg = mpc_modified.gen(idx, 2);
    Pmin = mpc_modified.gen(idx, 10);
    
    new_Pg = original_Pg - power_to_shed_per_unit;
    
    % 确保不会低于机组最小出力限制
    if new_Pg < Pmin
        new_Pg = Pmin;
    end
    
    mpc_modified.gen(idx, 2) = new_Pg;
    reduced_count = reduced_count + (original_Pg - new_Pg);
end

fprintf('信息: 常规发电机组总出力下调了 %.2f MW 来平衡新能源注入。\n', reduced_count);

end