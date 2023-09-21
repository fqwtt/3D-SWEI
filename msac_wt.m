function [bestModelParams] = msac_wt(allPoints, fitFunc, evalFunc, sampleSize, maxDistance)
    confidence = 99;
    numTrials  = 1000;
    numPts     = size(allPoints,1);
    idxTrial   = 1;
    bestModelParams  = evalin('base','model');
    [~, bestDis]    = evaluateModel(evalFunc, bestModelParams, allPoints, maxDistance);
    
    while idxTrial <= numTrials
        % 随机取点
        indices = randperm(numPts, sampleSize);
        samplePoints = allPoints(indices, :);
        modelParams = fitFunc(samplePoints);

        %计算loss
        [dis, accDis] = evaluateModel(evalFunc, modelParams, allPoints, maxDistance);

        % 更新参数
        if accDis < bestDis
            bestDis = accDis;
            bestModelParams = modelParams;
            inlierNum = sum(dis < maxDistance);
            num = computeLoopNumber(sampleSize, confidence, numPts, inlierNum);
            numTrials = min(numTrials, num);
        end
        idxTrial = idxTrial + 1;
    end
    
    % 用所有内点重新估计模型
    dis = evaluateModel(evalFunc, bestModelParams, allPoints, maxDistance);
    bestInliers = dis < maxDistance;
    bestModelParams = polyfit_wt(allPoints(bestInliers,1),allPoints(bestInliers,2)); 
end

%--------------------------------------------------------------------------
function [dis, accDis] = evaluateModel(evalFunc, modelIn, allPoints, maxDistance)
    dis = evalFunc(modelIn, allPoints);
    dis(dis > maxDistance) = maxDistance;
    accDis = sum(dis);
end

%--------------------------------------------------------------------------
function N = computeLoopNumber(sampleSize, confidence, pointNum, inlierNum)
    inlierProbability = (inlierNum/pointNum)^sampleSize;
    conf = 0.01 * confidence;
    num  = log10(1 - conf);
    den  = log10(1 - inlierProbability);
    N    = int32(ceil(num/den));
end

%--------------------------------------------------------------------------
function p = polyfit_wt(x,y)
    x = x(:);
    y = y(:);
    V(:,2) = ones(length(x),1);
    V(:,1) = x.*V(:,2);
    [Q,R] = qr(V);
    p = R\(Q'*y);
end
