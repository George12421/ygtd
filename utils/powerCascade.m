function success = powerCascade(mpc, initialFailureIndices, saveLocation, enableSaveData, flowType)

% --- 新增代码：计算并保存“第0阶段”（无故障）的状态 ---
% 1. 在移除任何线路之前，首先运行一次潮流计算
[mpc_stage0, success_stage0] = shedLoadAndSolvePF(mpc, flowType);
if(success_stage0 ~= 1)
    % 如果“第0阶段”都无法收敛，则直接返回失败
    success = success_stage0;
    return; 
end

% 2. 将“第0阶段”的状态作为 dataValues 的第一列
dataValues = appendSaveValues(mpc_stage0, []); 
% --- 新增代码结束 ---


% --- 原有代码开始（基于您的版本） ---
isOverLoaded = 1; 
% dataValues = []; % <-- 这行被移到新增代码中了
currentFailureIndices = initialFailureIndices;

while(isOverLoaded)
    % perform line failures (initial-failure or capacity-failure)
    % "第1阶段" 将在这里移除 initialFailureIndices
    mpc = removeInvalidLines(mpc, currentFailureIndices);
    
    % load shedding in each island and find power flow values
    [mpc, success] = shedLoadAndSolvePF(mpc, flowType);
    if(success ~= 1); return; end
    
    % record and save power values (追加“第1阶段”及以后的数据)
    dataValues = appendSaveValues(mpc, dataValues);
    
    % update the failure indices for next time step
    % [currentFailureIndices, isOverLoaded] = getOverLoadedLines(mpc);
    
    % --- 您原来的修改后的代码 ---
    
    % update the failure indices for next time step
    [allOverLoadedLines, isOverLoaded] = getOverLoadedLines(mpc);
    
    if (isOverLoaded)
        % --- 修改逻辑开始 ---
        % 原始代码将所有过载线路列表赋给 currentFailureIndices:
        % currentFailureIndices = allOverLoadedLines;
        
        % 新逻辑：为了实现“每一步只断一条线”，
        % 我们只从所有过载线路中选择 *第一条* 作为下个阶段的故障。
        % 其他过载线路（如果在新状态下仍然过载）将在后续的步骤中被检测到。
        
        currentFailureIndices = allOverLoadedLines(1);
        
        % 确保 isOverLoaded 仍然为 1，以继续 while 循环
        isOverLoaded = 1; 
        % --- 修改逻辑结束 ---
    else
        % 如果没有过载 (isOverLoaded = 0)，则正常结束
        currentFailureIndices = [];
    end
end

if(enableSaveData == 1)
    saveData(saveLocation, dataValues)
end

end