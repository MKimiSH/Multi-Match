function [affines, scores, matchConfigs] = K_MultiMatchDBSCAN(img_in, tpl_in, params)
% Use DBSCAN to first generate several good _configs_ clusters and then
% process them one by one (expand--evaluate--expand--...).
% Input: I, input image; T, template for matching
% Output: AA, array of Affine transformations; scores, corresponding match score for each affine
% 1. Downsample _I_ and _T_ w.r.t. some rules (same scale or to a proper scale and then adjust the scale range).
% 2. Use BBS to approximately locate possible matches. Generate a mask _M_, and only search where _M_==1 in the following steps.
% 2*. Run multiple BBS rounds with different scaling of _T_ may help constraint scale range.
% 3. Generate configs for searching.
% 4. Search for the 1st match _A1_.
% 4*. If _M_ consists of disjoint regions, then search each region for a best match.
% 5. Set _M(A1(T))=0_ and search for the 2nd match, and so on; stop when score is below a criterion.
% 6. Output the result.
% 暂时不考虑mask！
%% 0. Init
[BBSParams, matchParams] = CheckAllParams(params);
maxMatchRounds = 10;
useBigImages = 0;

%% 1. Downsample _I_ and _T_ w.r.t. some rules
[img, tpl, IResize, TResize] = AdjustSize(img_in, tpl_in, BBSParams.pz);
resizeFactor = IResize/TResize;
matchParams.resizeFactor = IResize/TResize;
matchParams.useBigImages = useBigImages;

%% 2. Use BBS to approximately locate possible matches. Generate a mask _M_
[M] = BBSPreprocess(img, tpl, BBSParams);

if useBigImages
    img = MakeOdd(img_in);
    tpl = MakeOdd(tpl_in);
    szI = size(img);
    szT = size(tpl);
    M = imresize(M, [szI(1), szI(2)]);
else
    szI = size(img);
    szT = size(tpl);
end

%% 2.5 Use SURF features to further restrain the rotation and scale search range
% the high-resolution version of the images are needed to get more features
[restrictedSearchRange, resM, scaleRotRange] = SearchRangeFromFeatureVoting(img_in, tpl_in, M); 
matchParams.searchRangeByFeature = restrictedSearchRange;
if ~isempty(scaleRotRange)
    if useBigImages
        matchParams.searchRange.minScale = scaleRotRange(:,1);
        matchParams.searchRange.maxScale = scaleRotRange(:,2);
    else
        matchParams.searchRange.minScale = scaleRotRange(:,1) * matchParams.resizeFactor;
        matchParams.searchRange.maxScale = scaleRotRange(:,2) * matchParams.resizeFactor;
    end
    matchParams.searchRange.minRotation = scaleRotRange(:,3);
    matchParams.searchRange.maxRotation = scaleRotRange(:,4);
end

%% 3. Generate configs for searching.
% M = ones(size(M)); % debug
[initConfigs, gridSize, matchParams] = InitializeConfigs(M, matchParams, szI, szT);
% 可以先initialize所有configs，但先根据feature得到的结果来搜索，用搜索结果去mask这些initConfigs。

%% 4. Search for matches guided by features
curAccept = 1;
matchRound = 1;
affines = [];
scores = zeros(maxMatchRounds, 1);
matchConfigs = zeros(maxMatchRounds, 6);
% matchRound = matchRound + 1;
A = [];
deltaFact = 1.511;
optMat = [];
% restrictedSearchRange = [];
for i = 1:size(restrictedSearchRange,1)
    
    [M, initConfigs] = MaskNewMatch(M, tpl, initConfigs, A);
    
    featSR = restrictedSearchRange(i,:);
    if ~useBigImages
        featSR(1:2) = featSR(1:2) * resizeFactor;
        featSR(5:8) = featSR(5:8) * IResize;
    end
    sr.minScale = featSR(1); % * resizeFactor;
    sr.maxScale = featSR(2); % * resizeFactor;
    sr.minRotation = featSR(3);
    sr.maxRotation = featSR(4);
    sr.minTx = featSR(5); % * IResize;
    sr.maxTx = featSR(6); % * IResize;
    sr.minTy = featSR(7); % * IResize;
    sr.maxTy = featSR(8); % * IResize;
    curMatchParams = matchParams;
    curMatchParams.searchRange = sr;
    
    [curInitConfigs, curGridSize, curMatchParams] = InitializeConfigs(M, curMatchParams, szI, szT);
    
    [A, score, config, accept] = FindOneMatch(img, tpl, M, curInitConfigs, curMatchParams, scores); % scores can be used as threshold
    
    curAccept = accept;
    config
    fprintf('Match #%d score (distance) = %.4f\n\n', matchRound, score);
    if(~curAccept)
        fprintf('Stop at match #%d!\n', matchRound);
    else
        affines = [affines, A];
        scores(matchRound) = score;
        matchConfigs(matchRound, :) = config;
        
        [optError,fullError,overlapError] = MatchingResult(tpl,img,A.tdata.T',optMat,sprintf('example %d', matchRound));
        fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
            matchRound, optError,256*optError,fullError,256*fullError,100*overlapError);
        matchRound = matchRound + 1;
    end
    
end

%% 5. Set _M(A1(T))=0_ and search for the 2nd match, and so on.
[goodConfigs, configClassList, numClasses] = FindGoodInitConfigs(img, tpl, M, initConfigs, matchParams);

while(curAccept && matchRound<=maxMatchRounds)
    fprintf('Match #%d, %d goodConfigs in this class: \n', matchRound, nnz(configClassList == matchRound));
    [M, initConfigs] = MaskNewMatch(M, tpl, initConfigs, A);
    initGoodConfigs = goodConfigs(configClassList == matchRound, :);
    curConfigs = ExpandConfigsRandom(initGoodConfigs, matchParams.steps, 1, 80, deltaFact);
    curConfigs = BoundConfigsTrSc(curConfigs, matchParams.bounds);
    [A, score, config, accept] = FindOneMatch(img, tpl, M, curConfigs, matchParams, scores); % scores can be used as threshold
    curAccept = accept;
    config
    fprintf('Match #%d score (distance) = %.4f\n\n', matchRound, score);
    if(~curAccept)
        fprintf('Stop at match #%d!\n', matchRound);
    else
        affines = [affines, A];
        scores(matchRound) = score;
        matchConfigs(matchRound, :) = config;
        
        [optError,fullError,overlapError] = MatchingResult(tpl,img,A.tdata.T',optMat,sprintf('example %d', matchRound));
        fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
            matchRound, optError,256*optError,fullError,256*fullError,100*overlapError);
    end
    matchRound = matchRound + 1;
end

%% 6. Output the result.
% 似乎不用做什么。
% 不对，得把大小还原回去……
affines = RecoverResizedAffines(affines, matchConfigs, IResize, TResize);

end

function [bbsparams, matchparams] = CheckAllParams(params)
% BBSParams: gamma, gz, rect(?)
% matchParams:
% epsilon, delta, photometricInvariance;
% searchRange: minTx, maxTx, minTy, maxTy, minRotation, maxRotation, minScale, maxScale

bbsparams = [];
matchparams = [];
if isfield(params, 'BBSParams'); bbsparams = params.BBSParams; end
if isfield(params, 'matchParams'); matchparams = params.matchParams; end


if(~isfield(bbsparams, 'gamma')); bbsparams.gamma = 2; end
if(~isfield(bbsparams, 'pz'));    bbsparams.pz = 3; end

if(~isfield(matchparams, 'epsilon')); matchparams.epsilon = 0.15; end
if(~isfield(matchparams, 'delta')); matchparams.delta = 0.25; end
if(~isfield(matchparams, 'photometricInvariance')); matchparams.photometricInvariance = 0; end
% searchRange will be initialized in InitializeConfigs()
if(~isfield(matchparams, 'searchRange')); matchparams.searchRange = []; end


end

function [affs] = RecoverResizedAffines(affs, configs, imgresize, tplresize)

configs(:, 1:2) = configs(:, 1:2) / imgresize; % tx, ty
configs(:, 4:5) = configs(:, 4:5) / (imgresize/tplresize); % sx, sy

for i = 1:length(affs)
%     aff = affs(i);
    
    affmat = CreateAffineTransformation(configs(i,:));
    affs(i) = maketform('affine', affmat');
end

end
