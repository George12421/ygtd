function initialFailureIndicesAll = getInitialFailures(numSamples, numLines, numFailures)
    % GETINITIALFAILURES: 按顺序且不重复地生成初始故障
    %
    % 描述:
    %   为 numSamples 个样本生成顺序且唯一的初始故障。
    %   例如：样本1 -> 故障线路1
    %         样本2 -> 故障线路2
    %         样本3 -> 故障线路3
    %
    %   Inputs:
    %       numSamples: 样本总数 (必须 <= numLines)
    %       numLines: 线路总数
    %       numFailures: 每次故障的线路数 (此版本只支持 1)
    %
    %   Outputs:
    %       initialFailureIndicesAll: 矩阵 (numFailures x numSamples)
    
    % --- 检查 1：此逻辑只为 N-1 (numFailures=1) 设计 ---
    if numFailures ~= 1
        error('“顺序且不重复”的逻辑目前只支持 numFailures = 1 的情况。');
    end
    
    % --- 检查 2：不重复样本数不能超过总线路数 ---
    if numSamples > numLines
        % 这是一个逻辑错误，我们必须阻止它
        error(['[getInitialFailures] 无法生成 %d 个“不重复”的故障样本，\n' ...
               '因为总共只有 %d 条线路。\n' ...
               '请将 numSamples 减少到 %d (或更低)，或者改用“顺序可重复”的逻辑。'], ...
               numSamples, numLines, numLines);
    end
    
    % --- 核心实现 (非常简单) ---
    
    % 1. 生成一个从 1 到 numSamples 的行向量
    %    例如 numSamples = 5 -> [1, 2, 3, 4, 5]
    %    由于上面的检查，我们保证这个序列中的每个数字都 <= numLines
    %    并且它们自动是“顺序且不重复的”。
    lineIndices = 1:numSamples;
    
    % 2. 直接赋值给输出矩阵
    initialFailureIndicesAll = zeros(numFailures, numSamples);
    initialFailureIndicesAll(1, :) = lineIndices;

end

% --- 注意：
% 原始版本中的子函数 getNewUniqueInitialFailure 不再需要了。