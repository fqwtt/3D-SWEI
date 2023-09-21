%RANSAC Fit a model to noisy data. 
%   [model, inlierIdx] = RANSAC(data, fitFcn, distFcn, sampleSize, maxDistance) 
%   fits a model to noisy data using the M-estimator SAmple Consensus 
%   (MSAC) algorithm, a version of RAndom SAmple Consensus algorithm.
%   
%   Inputs      Description 
%   ------      -----------
%   data        An M-by-N matrix, whose rows are data points to be modeled.
%               For example, for fitting a line to 2-D points, data would be an 
%               M-by-2 matrix of [x,y] coordinates. For fitting a geometric 
%               transformation between two sets of matched 2-D points, 
%               the coordinates can be concatenated into an M-by-4 matrix.
%   
%   fitFcn      A handle to a function, which fits the model to a minimal
%               subset of data. The function must be of the form 
%                 model = fitFcn(data)
%               model returned by fitFcn can be a cell array, if it is
%               possible to fit multiple models to the data. 
%  
%   distFcn     A handle to a function, which computes the distances from the
%               model to the data. The function must be of the form
%                 distances = distFcn(model, data)
%               If model is an N-element cell array, then distances must be an 
%               M-by-N matrix. Otherwise, distances must be an M-by-1 vector.
%  
%   sampleSize  A positive integer scalar containing the minimum size of a 
%               sample from data required by fitFcn to fit a model.
%
%   maxDistance A positive numeric scalar specifying the distance threshold 
%               for finding outliers. Increasing this value will make the 
%               algorithm converge faster, but may adversely affect the 
%               accuracy of the result.
%  
%   Outputs     Description
%   -------     -----------
%   model       The model, which best fits the data.
% 
%   inlierIdx   An M-by-1 logical array, specifying which data points are 
%               inliers.
%
%   [..., status] = RANSAC(...) additionally returns a status code. If the 
%   status output is not specified, the function will issue an error if 
%   the number of data points is less than sampleSize, or if a model cannot
%   be estimated. The status can have the following values:
%      0: No error.
%      1: The number of data points is less than sampleSize.
%      2: The number of inliers is less than sampleSize.
%  
%   [...] = RANSAC(..., Name, Value) specifies additional name-value pair
%   arguments described below:
%  
%   'ValidateModelFcn'    Handle to a function, which returns true if the 
%                         model is valid and false otherwise. Certain subsets
%                         of data may be degenerate, causing fitFcn to 
%                         return an invalid model. The function must be of
%                         the form,
%                                  isValid = validateModelFcn(model)
%
%                         If ValidateModelFcn is not specified, all models
%                         returned by fitFcn are assumed to be valid.
% 
%   'MaxSamplingAttempts' Positive integer scalar specifying the maximum
%                         number of attempts to find a sample, which yields
%                         a valid model. This parameters is used only if 
%                         ValidateModelFun is set.
%  
%                         Default: 100
%  
%   'MaxNumTrials'        Positive integer scalar specifying the maximum 
%                         number of random trials for finding the best model.
%                         The actual number of trials depends on the data, 
%                         and the values of the MaxDistance and Confidence 
%                         parameters. Increasing this value will improve 
%                         robustness of the output at the expense of 
%                         additional computation.
%  
%                         Default: 1000
%  
%   'Confidence'          Scalar value greater than 0 and less than 100.
%                         Specifies the desired confidence (in percentage)
%                         for finding the maximum number of inliers.
%                         Increasing this value will improve the robustness
%                         of the output at the expense of additional
%                         computation.
%  
%                         Default: 99
%
%  Class Support
%  -------------
%  data can be double or single. fitFcn and distFcn must be function
%  handles. sampleSize is a numeric scalar. maxDistance can be double or
%  single.
%
%  Example: Fit a line to set of 2-D points.
%  -----------------------------------------
%  % Load and plot a set of noisy 2D points.
%  load 'pointsForLineFitting.mat';
%  plot(points(:,1), points(:,2), '*');
%  hold on
%
%  % Fit a line using linear least squares.
%  modelLeastSquares = polyfit(points(:,1), points(:,2), 1);
%  x = [min(points(:,1)), max(points(:,1))];
%  y = modelLeastSquares(1)*x + modelLeastSquares(2);
%  plot(x, y, 'r-');
%
%  % Fit a line to the points using M-estimator SAmple Consensus algorithm.
%  sampleSize = 2;
%  maxDistance = 2;
%  fitLineFcn  = @(points) polyfit(points(:,1), points(:,2), 1);
%  evalLineFcn = ...
%     @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2, 2);
%     
%  [modelRANSAC, inlierIdx] = RANSAC(points, fitLineFcn, evalLineFcn, ...
%     sampleSize, maxDistance);
%
%  % Re-fit a line to the inliers.
%  modelInliers = polyfit(points(inlierIdx,1), points(inlierIdx,2), 1);
%
%  % Display the line.
%  inlierPts = points(inlierIdx, :);
%  x = [min(inlierPts(:,1)); max(inlierPts(:,1))];
%  y = modelInliers(1)*x + modelInliers(2);
%  plot(x, y, 'g-');
%  legend('Noisy Points', 'Least Squares Fit', 'Robust Fit');
%  hold off
%
%  See also estimateEssentialMatrix.

% References:
%   P. H. S. Torr and A. Zisserman, "MLESAC: A New Robust Estimator
%   with Application to Estimating Image Geometry," Computer Vision
%   and Image Understanding, 2000.

% Copyright 2016-2021 The MathWorks, Inc.

%#codegen

function [model, inlierIdx, status] = ransac1(data, fitFun, distFun, sampleSize,...
    maxDistance, varargin)

if isSimMode
    if nargin > 5
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
end

[params, funcs] = parseInputs(data, fitFun, distFun, sampleSize, ...
    maxDistance, varargin{:});

% List of status codes
statusCode = struct(...
    'NoError',           int32(0),...
    'NotEnoughPts',      int32(1),...
    'NotEnoughInliers',  int32(2));

reachedMaxSkipTrials = false;

if size(data, 1) < sampleSize
    status = statusCode.NotEnoughPts;
    model = cast([], 'like', data);
    inlierIdx = false(size(data, 1), 1);
else    
    [isFound, model, inlierIdx, reachedMaxSkipTrials] = ...
        msac1(data, params, funcs);
        
    if isFound
        status = statusCode.NoError;
    else
        status = statusCode.NotEnoughInliers;
    end
end

if reachedMaxSkipTrials
    coder.internal.warning('vision:ransac:reachedMaxSkipTrials', ...
        params.maxSkipTrials);
end

if nargout < 3
    checkRuntimeStatus(statusCode, status, sampleSize);
end

%==========================================================================
% Check runtime status and report error if there is one
%==========================================================================
function checkRuntimeStatus(statusCode, status, sampleSize)
coder.internal.errorIf(status==statusCode.NotEnoughPts, ...
    'vision:ransac:notEnoughDataPts', sampleSize);

coder.internal.errorIf(status==statusCode.NotEnoughInliers, ...
    'vision:ransac:notEnoughInliers');

%--------------------------------------------------------------------------
function [params, funcs] = parseInputs(data, fitFun, distFun,...
    sampleSize, maxDistance, varargin)

coder.internal.prefer_const( varargin{:} );

validateattributes(data, {'single', 'double'}, ...
    {'nonempty', '2d', 'real', 'nonsparse'}, mfilename, 'data');

validateattributes(fitFun,  {'function_handle'}, {'scalar'}, mfilename, 'fitFun');
validateattributes(distFun, {'function_handle'}, {'scalar'}, mfilename, 'distFun');

validateattributes(sampleSize, {'numeric'}, {'scalar', 'integer', 'positive'}, ...
    mfilename, 'sampleSize');
checkMaxDistance(maxDistance);

funcs.fitFunc   = fitFun;
funcs.evalFunc  = distFun;
params.sampleSize  = sampleSize;
params.maxDistance = maxDistance;
params.recomputeModelFromInliers = false;

defaults = struct('MaxNumTrials', {1000}, 'Confidence', {99}, ...
    'ValidateModelFcn', {@validateModelDefault}, 'FitFcnParameters', {{}}, ...
    'MaxSamplingAttempts', {100});

if isSimMode
    parser = inputParser;
    parser.CaseSensitive = false;
    parser.FunctionName = mfilename;
    parser.addParameter('ValidateModelFcn',    defaults.ValidateModelFcn,    @checkValidateModelFcn);
    parser.addParameter('MaxNumTrials',        defaults.MaxNumTrials,        @checkMaxNumTrials);
    parser.addParameter('MaxSamplingAttempts', defaults.MaxSamplingAttempts, @checkMaxSamplingAttempts);
    parser.addParameter('Confidence',          defaults.Confidence,          @checkConfidence);
    
    parser.parse(varargin{:});
    params.confidence    = parser.Results.Confidence;
    params.maxNumTrials  = parser.Results.MaxNumTrials;
    params.maxSkipTrials = parser.Results.MaxSamplingAttempts;
    
    funcs.checkFunc = parser.Results.ValidateModelFcn;


else
    % Define parser mapping struct
    pvPairs = struct( ...
        'MaxNumTrials',       uint32(0), ...
        'ValidateModelFcn',   false, ...
        'FitFcnParameters',   uint32(0), ...
        'MaxSamplingAttempts',uint32(0), ...
        'Confidence',         uint32(0));
    % Specify parser options
    poptions = struct( ...
        'CaseSensitivity', false, ...
        'StructExpand',    true, ...
        'PartialMatching', true);
    
    % Parse PV pairs
    pstruct = coder.internal.parseParameterInputs(pvPairs, ...
        poptions, varargin{:});
    % Extract inputs
    funcs.checkFunc         = coder.internal.getParameterValue(pstruct.ValidateModelFcn, defaults.ValidateModelFcn, varargin{:});
    params.confidence       = coder.internal.getParameterValue(pstruct.Confidence, defaults.Confidence, varargin{:});
    params.maxNumTrials     = coder.internal.getParameterValue(pstruct.MaxNumTrials, defaults.MaxNumTrials, varargin{:});
    params.maxSkipTrials    = coder.internal.getParameterValue(pstruct.MaxSamplingAttempts, defaults.MaxSamplingAttempts, varargin{:});
    
    % Validate inputs
    checkValidateModelFcn(funcs.checkFunc);
    checkMaxNumTrials(params.maxNumTrials);
    checkMaxSamplingAttempts(params.maxSkipTrials);
    checkConfidence(params.confidence);
end

%--------------------------------------------------------------------------
function tf = validateModelDefault(varargin)
tf = true;

%--------------------------------------------------------------------------
function checkValidateModelFcn(fun)
validateattributes(fun, {'function_handle'}, {'scalar'}, mfilename, 'ValidateModelFcn');

%--------------------------------------------------------------------------
function checkMaxNumTrials(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'integer', 'positive'}, mfilename, ...
    'MaxNumTrials');

%--------------------------------------------------------------------------
function checkMaxSamplingAttempts(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'integer', 'positive'}, mfilename, ...
    'MaxSamplingAttempts');

%--------------------------------------------------------------------------
function checkConfidence(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'positive', '<', 100}, mfilename, ...
    'Confidence');

%--------------------------------------------------------------------------
function checkMaxDistance(value)
validateattributes(value,{'single','double'}, ...
    {'real', 'nonsparse', 'scalar','nonnegative','finite'}, mfilename, ...
    'maxDistance');

%--------------------------------------------------------------------------
function flag = isSimMode()

flag = isempty(coder.target);


%--------------------------------------------------------------------------
function [isFound, bestModelParams, inliers, reachedMaxSkipTrials] = msac1(...
    allPoints, params, funcs, varargin)
% MSAC M-estimator SAmple Consensus (MSAC) algorithm that is used for point
% cloud model fitting. allPoints must be an M-by-N matrix, where each point
% is a row vector.
%
% allPoints - M-by-2 or M-by-2-by-2 array of [x y] coordinates
%
% params    - struct containing the following fields:
%               sampleSize
%               maxDistance
%               confidence
%               maxNumTrials
%
% funcs     - struct containing the following function handles
%               fitFunc
%               evalFunc
%               checkFunc

% Copyright 2015-2020 The MathWorks, Inc.
%
% References:
% ----------
%   P. H. S. Torr and A. Zisserman, "MLESAC: A New Robust Estimator with
%   Application to Estimating Image Geometry," Computer Vision and Image
%   Understanding, 2000.

%#codegen

confidence = params.confidence;
sampleSize = params.sampleSize;
maxDistance = params.maxDistance;

threshold = cast(maxDistance, 'like', allPoints);
numPts    = size(allPoints,1);
idxTrial  = 1;
numTrials = int32(params.maxNumTrials);
maxDis    = cast(threshold * numPts, 'like', allPoints);
bestDis   = maxDis;

if isfield(params, 'defaultModel')
    bestModelParams = params.defaultModel;
else
    bestModelParams = zeros(0, 'like', allPoints);
end

if isfield(params, 'maxSkipTrials')
    maxSkipTrials = params.maxSkipTrials;
else
    maxSkipTrials = params.maxNumTrials * 10;
end
skipTrials = 0;
reachedMaxSkipTrials = false;

bestInliers = false(numPts, 1);

% Create a random stream. It uses a fixed seed for the testing mode and a
% random seed for other mode.
coder.extrinsic('vision.internal.testEstimateGeometricTransform');
if isempty(coder.target) && vision.internal.testEstimateGeometricTransform
    rng('default');
end


while idxTrial <= numTrials && skipTrials < maxSkipTrials
    % Random selection without replacement
    indices = randperm(numPts, sampleSize);
    
    % Compute a model from samples
    samplePoints = allPoints(indices, :, :);
    modelParams = funcs.fitFunc(samplePoints, varargin{:});
    
    % Validate the model
    isValidModel = funcs.checkFunc(modelParams, varargin{:});
    
    if isValidModel
        % Evaluate model with truncated loss
        [model, dis, accDis] = evaluateModel(funcs.evalFunc, modelParams, ...
            allPoints, threshold, varargin{:});
        
        % Update the best model found so far
        if accDis < bestDis
            bestDis = accDis;
            bestInliers = dis < threshold;
            bestModelParams = model;
            inlierNum = cast(sum(dis < threshold), 'like', allPoints);
%             num = vision.internal.ransac.computeLoopNumber(sampleSize, ...
%                 confidence, numPts, inlierNum);
%             numTrials = min(numTrials, num);
        end
        
        idxTrial = idxTrial + 1;
    else
        skipTrials = skipTrials + 1;
    end
end

isFound = funcs.checkFunc(bestModelParams(:), varargin{:}) && ...
    ~isempty(bestInliers) && sum(bestInliers(:)) >= sampleSize;
if isFound
    if isfield(params, 'recomputeModelFromInliers') && ...
            params.recomputeModelFromInliers
        modelParams = funcs.fitFunc(allPoints(bestInliers, :, :), varargin{:});
        [bestModelParams, dis] = evaluateModel(funcs.evalFunc, modelParams, ...
            allPoints, threshold, varargin{:});
        isValidModel = funcs.checkFunc(bestModelParams(:), varargin{:});
        inliers = (dis < threshold);
        if ~isValidModel || ~any(inliers)
            isFound = false;
            inliers = false(size(allPoints, 1), 1);
            return;
        end
    else
        inliers = bestInliers;
    end
    
%     if numTrials >= int32(params.maxNumTrials)
%         coder.internal.warning('vision:ransac:maxTrialsReached');
%     end
else
    inliers = false(size(allPoints, 1), 1);
end

reachedMaxSkipTrials = skipTrials >= maxSkipTrials;

%--------------------------------------------------------------------------
function [modelOut, distances, sumDistances] = evaluateModel(evalFunc, modelIn, ...
    allPoints, threshold, varargin)
dis = evalFunc(modelIn, allPoints, varargin{:});
dis(dis > threshold) = threshold;
accDis = sum(dis);
if iscell(modelIn)
    [sumDistances, minIdx] = min(accDis);
    distances = dis(:, minIdx);
    modelOut = modelIn{minIdx(1)};
else
    distances = dis;
    modelOut = modelIn;
    sumDistances = accDis;
end


