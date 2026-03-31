function initialFailureIndicesAll = getInitialFailures_v1(numSamples, numLines, numFailures)
    % GETINITIALFAILURES: Gets a list of unique initial failure contingencies
    %
    %   Inputs:
    %       numSamples: Number of samples to generate
    %       numLines: Total number of lines in the grid
    %       numFailures: Number of initial failures (1 or 2)
    %
    %   Outputs:
    %       initialFailureIndicesAll: A matrix (numFailures x numSamples)
    
    initialFailureIndicesAll = zeros(numFailures, numSamples);
    
    % uniquenessTracker 确保我们不会生成重复的故障组合
    % 我们将使用一个 2D 矩阵。
    % - 对于 N-1 故障 (numFailures = 1)，我们将使用 (index, 1) 来标记。
    % - 对于 N-2 故障 (numFailures = 2)，我们将使用 (index1, index2) 来标记。
    uniquenessTracker = zeros(numLines, numLines);
    
    for jj = 1:numSamples
        % 调用修复后的子函数
        [failureIndices, uniquenessTracker] = getNewUniqueInitialFailure(numLines, numFailures, uniquenessTracker);
        initialFailureIndicesAll(:,jj) = failureIndices;
    end
end

% --- 这是包含核心修复的子函数 ---
function [failureIndices, uniquenessTracker] = getNewUniqueInitialFailure(numLines, numFailures, uniquenessTracker)
    
    % 生成一个随机的故障组合
    failureIndices = randperm(numLines, numFailures);
    failureIndices = sort(failureIndices);
    
    % +++ 关键修复：检查 numFailures +++
    if numFailures == 1
        % --- N-1 (单一故障) 的处理逻辑 ---
        
        % 检查这个单一故障是否已经被选过
        % (我们利用 uniquenessTracker 的第一列来做这件事)
        while(uniquenessTracker(failureIndices(1), 1) == 1)
            failureIndices = randperm(numLines, numFailures);
            % (numFailures=1 时不需要排序)
        end
        % 标记这个故障已被使用
        uniquenessTracker(failureIndices(1), 1) = 1;
        
    elseif numFailures == 2
        % --- N-2 (双重故障) 的原始逻辑 ---
        
        % 检查这对故障组合是否已经被选过
        while(uniquenessTracker(failureIndices(1), failureIndices(2)) == 1)
            failureIndices = randperm(numLines, numFailures);
            failureIndices = sort(failureIndices);
        end
        % 标记这对故障已被使用
        uniquenessTracker(failureIndices(1), failureIndices(2)) = 1;
        uniquenessTracker(failureIndices(2), failureIndices(1)) = 1; % 确保对称
        
    else
        % --- 对其他情况抛出错误 ---
        error('getInitialFailures.m 只被设计用于处理 1 或 2 个初始故障。');
    end
    
    failureIndices = failureIndices'; % 转换为列向量
end