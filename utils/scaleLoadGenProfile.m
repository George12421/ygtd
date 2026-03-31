function mpc = scaleLoadGenProfile(mpc, loadScale)
    % ！！删除这里的 scaleLimitLow 和 scaleLimitHigh ！！

    loadVals = mpc.bus(:,3);
    genVals = mpc.gen(:,2);

    if(strcmp(loadScale, 'random'))
        % --- 将定义移动到 IF 内部 ---
        scaleLimitLow = mpc.scaleLimitLow;
        scaleLimitHigh = mpc.scaleLimitHigh;
        % ---
        loadScale = scaleLimitLow + (scaleLimitHigh - scaleLimitLow)*(rand());
    %     loadScale = randsample([1.00:0.1:2.00], 1);
    end

    if(strcmp(loadScale, 'random-nonuniform'))
        % --- 将定义移动到 IF 内部 ---
        scaleLimitLow = mpc.scaleLimitLow;
        scaleLimitHigh = mpc.scaleLimitHigh;
        % ---
        
        % (注意: 这里的 'if(rand() < 0)' 条件永远不会成立)
        if(rand() < 0) 
            loadScale = scaleLimitLow +(scaleLimitHigh - scaleLimitLow)*(rand());
        else
            loadScale = scaleLimitLow +(scaleLimitHigh - scaleLimitLow)*(rand(size(loadVals)));
        end
    end

    if(strcmp(loadScale, 'non-uniform'))
        loadScale = mpc.knownScaling;
    end

    % 这部分对于 loadScale = 1.0 是完全正常的
    loadVals = loadScale.*loadVals;
    genVals  = sum(loadVals)/sum(genVals)*genVals;
    mpc.bus(:,3) = loadVals;
    mpc.gen(:,2) = genVals;
end