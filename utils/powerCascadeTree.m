function all_sequences = powerCascadeTree(mpc, initialFailureIndices, flowType, max_depth)
    % 生成树状级联序列的主入口
    % 返回值 all_sequences: 包含所有可能的分支序列的 Cell 数组
    all_sequences = {}; 
    
    % --- 第 0 阶段：初始无故障稳态 ---
    [mpc_stage0, success_stage0] = shedLoadAndSolvePF(mpc, flowType);
    if(success_stage0 ~= 1)
        return; % 初始状态不收敛，直接退出
    end
    dataValues_initial = appendSaveValues(mpc_stage0, []); 
    
    % --- 第 1 阶段：断开初始故障 ---
    mpc_stage1 = removeInvalidLines(mpc_stage0, initialFailureIndices);
    [mpc_stage1, success_stage1] = shedLoadAndSolvePF(mpc_stage1, flowType);
    if(success_stage1 ~= 1)
        return; 
    end
    dataValues_stage1 = appendSaveValues(mpc_stage1, dataValues_initial);
    
    % --- 开始递归探索树状级联 (当前深度为 1) ---
    all_sequences = dfs_explore(mpc_stage1, dataValues_stage1, flowType, 1, max_depth, all_sequences);
end

% ----------------- 内部递归函数 (用于探索所有分支) -----------------
function seqs = dfs_explore(mpc_current, dataValues_current, flowType, current_depth, max_depth, seqs)
    
    % 1. 如果达到了设定的最大级联深度 M，终止当前分支并保存
    if current_depth >= max_depth
        seqs{end+1} = dataValues_current; 
        return;
    end
    
    % 2. 检测当前状态下哪些线路过载了
    [overloadedLines, isOverLoaded] = getOverLoadedLines(mpc_current);
    
    % 3. 如果没有线路过载，该分支级联自然停止，保存序列
    if ~isOverLoaded
        seqs{end+1} = dataValues_current;
        return;
    end
    
    % 4. 【分支遍历 (Branching)】：对每一条过载线路单独移除，探索多条可能序列
    for i = 1:length(overloadedLines)
        line_to_remove = overloadedLines(i);
        
        % 复制当前电网状态，确保分支之间互不干扰
        mpc_next = mpc_current; 
        
        % 仅移除这一条过载线路
        mpc_next = removeInvalidLines(mpc_next, line_to_remove);
        [mpc_next, success] = shedLoadAndSolvePF(mpc_next, flowType);
        
        if success == 1
            % 记录这一步的数据，进入下一层递归
            dataValues_next = appendSaveValues(mpc_next, dataValues_current);
            seqs = dfs_explore(mpc_next, dataValues_next, flowType, current_depth + 1, max_depth, seqs);
        else
            % 如果导致潮流不收敛(系统崩溃)，该分支到此结束，保存崩塌前的状态
            seqs{end+1} = dataValues_current; 
        end
    end
end