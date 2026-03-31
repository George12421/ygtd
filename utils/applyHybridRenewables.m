% Filename: applyHybridRenewables.m
% Location: matlabmodels/utils/
% (这是“混合能源、多节点注入”的最终版)

function mpc_modified = applyHybridRenewables(mpc_original, solar_scenario, wind_scenario)
%% applyHybridRenewables: 将【光伏和风电】的出力场景应用到MATPOWER case中
%
%   Inputs:
%       mpc_original: 原始的MATPOWER case结构体
%       solar_scenario: 一个结构体，包含光伏的 .bus_numbers 和 .power_mw_list
%       wind_scenario:  一个结构体，包含风电的 .bus_numbers 和 .power_mw_list
%       conventional_gens: (我们将其移到solar_scenario结构体内部)
%
%   Outputs:
%       mpc_modified: 应用新能源场景后，修改过的mpc结构体

% 深度复制原始mpc对象
mpc_modified = mpc_original;
P_new_renewable_total = 0; % 总注入功率计数器

% --- 步骤 1: 注入所有【光伏】电力 ---
if isfield(solar_scenario, 'bus_numbers') && ~isempty(solar_scenario.bus_numbers)
    if length(solar_scenario.bus_numbers) ~= length(solar_scenario.power_mw_list)
        error('错误: 光伏母线列表和功率列表的长度必须一致！');
    end
    
    for i = 1:length(solar_scenario.bus_numbers)
        current_bus = solar_scenario.bus_numbers(i);
        current_power = solar_scenario.power_mw_list(i);
        
        bus_idx = find(mpc_modified.bus(:, 1) == current_bus);
        if isempty(bus_idx)
            error('错误: 在案例中未找到指定的光伏母线编号 %d', current_bus);
        end
        
        original_load_pd = mpc_modified.bus(bus_idx, 3);
        mpc_modified.bus(bus_idx, 3) = original_load_pd - current_power;
        P_new_renewable_total = P_new_renewable_total + current_power;
    end
    fprintf('信息: 注入 %.2f MW 光伏电力。\n', sum(solar_scenario.power_mw_list));
end

% --- 步骤 2: 注入所有【风电】电力 ---
if isfield(wind_scenario, 'bus_numbers') && ~isempty(wind_scenario.bus_numbers)
    if length(wind_scenario.bus_numbers) ~= length(wind_scenario.power_mw_list)
        error('错误: 风电母线列表和功率列表的长度必须一致！');
    end
    
    for i = 1:length(wind_scenario.bus_numbers)
        current_bus = wind_scenario.bus_numbers(i);
        current_power = wind_scenario.power_mw_list(i);
        
        bus_idx = find(mpc_modified.bus(:, 1) == current_bus);
        if isempty(bus_idx)
            error('错误: 在案例中未找到指定的风电母线编号 %d', current_bus);
        end
        
        original_load_pd = mpc_modified.bus(bus_idx, 3);
        mpc_modified.bus(bus_idx, 3) = original_load_pd - current_power;
        P_new_renewable_total = P_new_renewable_total + current_power;
    end
    fprintf('信息: 注入 %.2f MW 风力电力。\n', sum(wind_scenario.power_mw_list));
end

fprintf('信息: 总共注入 %.2f MW 新能源。\n', P_new_renewable_total);

% --- 步骤 3: 调整常规发电机出力以维持系统总功率平衡 ---
% (我们假设常规机组列表在 solar_scenario 中定义)
conventional_gens_indices = solar_scenario.conventional_gens; 

% (这部分平衡逻辑与之前完全相同)
if any(conventional_gens_indices > size(mpc_modified.gen, 1))
    error('错误: 常规发电机索引超出范围。');
end
total_conventional_power = sum(mpc_modified.gen(conventional_gens_indices, 2));
total_conventional_pmin = sum(mpc_modified.gen(conventional_gens_indices, 10));

if (total_conventional_power - P_new_renewable_total) < total_conventional_pmin
    fprintf('警告: 新能源总注入量过大，常规机组下调后将低于其总最小出力。\n');
    P_to_reduce = total_conventional_power - total_conventional_pmin;
else
    P_to_reduce = P_new_renewable_total;
end

power_to_shed_per_unit = P_to_reduce / length(conventional_gens_indices);
reduced_count = 0;
for i = 1:length(conventional_gens_indices)
    idx = conventional_gens_indices(i);
    original_Pg = mpc_modified.gen(idx, 2);
    Pmin = mpc_modified.gen(idx, 10);
    new_Pg = original_Pg - power_to_shed_per_unit;
    if new_Pg < Pmin; new_Pg = Pmin; end
    mpc_modified.gen(idx, 2) = new_Pg;
    reduced_count = reduced_count + (original_Pg - new_Pg);
end
fprintf('信息: 常规发电机组总出力下调了 %.2f MW 来平衡。\n', reduced_count);

end