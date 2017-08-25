function [affines, scores, matchConfigs, timeList, fullErrors, nFeatMatch, nOtherMatch] = K_MultiMatchFeaturePoints(img_in, tpl_in, params)
% Use *SURF* to first generate several good _configs_ clusters and then
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

useBigImages = 1; % use original images or subsampled images for matching (default=1)
useClustering = 0; % use DBSCAN clustering or not (default=0)

[BBSParams, matchParams] = CheckAllParams(params);
maxMatchRounds = params.nMatches;
useBBS = params.useBBS;
useSURF = params.useSURF;

timeList.BBSTime = 0;
timeList.featureTime = 0;
timeList.clusteringTime = 0;
timeList.featureMatchTime = [];
timeList.otherMatchTime = [];
fullErrors = [];

%% 1. Downsample _I_ and _T_ w.r.t. some rules
[img, tpl, IResize, TResize] = AdjustSize(img_in, tpl_in, BBSParams.pz);
resizeFactor = IResize/TResize;
matchParams.resizeFactor = IResize/TResize;
matchParams.useBigImages = useBigImages;

%% 2. Use BBS to approximately locate possible matches. Generate a mask _M_
if(useBBS)
    [M, timeList.BBSTime] = BBSPreprocess(img, tpl, BBSParams);
    figure, imshow(M);
else
    M = ones(size(img, 1), size(img, 2));
end

if useBigImages
    img = MakeOdd(img_in);
    tpl_in = imresize(tpl_in, 1/resizeFactor);
    tpl = MakeOdd(tpl_in);
    szI = size(img);
    szT = size(tpl);
    M = imresize(M, [szI(1), szI(2)]);
else
    szI = size(img);
    szT = size(tpl);
end

%% 3. Use SURF features to further restrain the rotation and scale search range
% the high-resolution version of the images are needed to get more features
if(useSURF)
    [restrictedSearchRange, resM, scaleRotRange, ~, ~, timeList.featureTime] = SearchRangeFromFeatureVoting(img_in, tpl_in, M);
    matchParams.searchRangeByFeature = restrictedSearchRange;
else
    restrictedSearchRange = [];
end

%% 4. Generate configs for searching.
% M = ones(size(M)); % debug
[initConfigs, gridSize, matchParams] = InitializeConfigs(M, matchParams, szI, szT);
% 可以先initialize所有configs，但先根据feature得到的结果来搜索，用搜索结果去mask这些initConfigs。

%% 5. Search for matches guided by features
curAccept = 1;
matchRound = 1;
nFeatMatch = 0;
nOtherMatch = 0;
affines = [];
scores = []; %zeros(maxMatchRounds, 1);
matchConfigs = zeros(maxMatchRounds, 6);
% matchRound = matchRound + 1;
A = [];
deltaFact = 1.511;
optMat = [];
% restrictedSearchRange = [];
scalerotRanges = [];
acceptFeatureRange = [];
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
    if isempty(curInitConfigs)
        fprintf('Feature based match #%d not found because of no configs!!\n', matchRound);
        acceptFeatureRange = [acceptFeatureRange, 0];
        continue;
    end
    [A, score, config, accept, curtime] = FindOneMatch(img, tpl, M, curInitConfigs, curMatchParams, scores); 
    timeList.featureMatchTime = [timeList.featureMatchTime, curtime];
    curAccept = accept;
    acceptFeatureRange = [acceptFeatureRange, accept];
    config
    fprintf('Match #%d score (distance) = %.4f\n\n', matchRound, score);
    if(~curAccept)
        fprintf('Feature based match #%d not found because of not accepted!!\n', matchRound);
    else
        scalerotRanges = [scalerotRanges; featSR];
        affines = [affines, A];
        scores = [scores, score];
        matchConfigs(matchRound, :) = config;
        
%         [optError,fullError,overlapError] = MatchingResult(tpl,img,A.tdata.T',optMat,sprintf('example %d', matchRound));
%         fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
%             matchRound, optError,256*optError,fullError,256*fullError,100*overlapError);
%         fullErrors = [fullErrors, fullError];
        matchRound = matchRound + 1;
    end
    
end
nFeatMatch = matchRound - 1;

% draw feature voting search ranges
% figure, imshow(img_in); hold on
% for j = 1:length(acceptFeatureRange)
%     if(acceptFeatureRange(j))
%         rectangle('Position', [restrictedSearchRange(j,5) + size(img_in,2)/2, ...
%             restrictedSearchRange(j,7) + size(img_in,1)/2, restrictedSearchRange(j,6) - restrictedSearchRange(j,5), ...
%             restrictedSearchRange(j,8) - restrictedSearchRange(j,7) ], 'EdgeColor', 'g', 'LineWidth', 1.5);
%     else
%         rectangle('Position', [restrictedSearchRange(j,5) + size(img_in,2)/2, ...
%             restrictedSearchRange(j,7) + size(img_in,1)/2, restrictedSearchRange(j,6) - restrictedSearchRange(j,5), ...
%             restrictedSearchRange(j,8) - restrictedSearchRange(j,7) ], 'EdgeColor', 'y', 'LineWidth', 1.5);
%     end
% end

if ~isempty(scalerotRanges)
    allmins = min(scalerotRanges(:,1));
    allmaxs = max(scalerotRanges(:,2));
    allminr = min(scalerotRanges(:,3));
    allmaxr = max(scalerotRanges(:,4));
    scaleRotRange = [allmins, allmaxs, allminr, allmaxr]
else
    scaleRotRange = [];
end

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


%% 6. Set _M(A1(T))=0_ and search for the remaining matches.
if(nFeatMatch >= maxMatchRounds)
    fprintf('\n!!!!!!!!ENDING because enough feature matches!!!!!!!!!!\n');
    fprintf('nFeatMatch = %d, nOtherMatch = 0\n', nFeatMatch);
    return;
end

[initConfigs, gridSize, matchParams] = InitializeConfigs(M, matchParams, szI, szT);

% assess initConfigs and find good ones in it to start further matching
[goodConfigs, configClassList, numClasses, allDistances, initConfigs, initTime] = FindGoodInitConfigs(img, tpl, M, initConfigs, matchParams, useClustering);
timeList.otherMatchTime = [timeList.otherMatchTime, initTime];
newMatchRound = 1;

if ~useClustering
    initGoodConfigs = goodConfigs(configClassList == 1, :);
    curConfigs = ExpandConfigsRandom(initGoodConfigs, matchParams.steps, 1, 80, deltaFact);
    curConfigs = BoundConfigsTrSc(curConfigs, matchParams.bounds);
end
curAccept = 1;
while(curAccept && matchRound<=maxMatchRounds)
    if useClustering
        fprintf('\nMatch #%d, %d goodConfigs in this class: \n', matchRound, nnz(configClassList == matchRound));
        [M, initConfigs] = MaskNewMatch(M, tpl, initConfigs, A);
        initGoodConfigs = goodConfigs(configClassList == newMatchRound, :);
        curConfigs = ExpandConfigsRandom(initGoodConfigs, matchParams.steps, 1, 80, deltaFact);
        curConfigs = BoundConfigsTrSc(curConfigs, matchParams.bounds);
        newMatchRound = newMatchRound + 1;
    else
        [M, curConfigs] = MaskNewMatch(M, tpl, curConfigs, A);
        fprintf('\n--->>>Match #%d\n', matchRound);
        if(size(curConfigs, 1) < 50) % 补充一下资源
            [initGoodConfigs, ~, ~, allDistances, initConfigs, initTime] = ...
                FindGoodInitConfigs(img, tpl, M, initConfigs, matchParams, useClustering, allDistances);
            curConfigs = ExpandConfigsRandom(initGoodConfigs, matchParams.steps, 1, 80, deltaFact);
            curConfigs = BoundConfigsTrSc(curConfigs, matchParams.bounds);
        end
    end
        
    [A, score, config, accept, curtime] = FindOneMatch(img, tpl, M, curConfigs, matchParams, scores); % scores can be used as threshold
    timeList.otherMatchTime = [timeList.otherMatchTime, curtime];
    curAccept = accept;
    config
    fprintf('Match #%d score (distance) = %.4f\n\n', matchRound, score);
    if(~curAccept)
        fprintf('Stop at match #%d!\n', matchRound);
        break;
    else
        affines = [affines, A];
        scores = [scores, score];
        matchConfigs(matchRound, :) = config;
        
%         [optError,fullError,overlapError] = MatchingResult(tpl,img,A.tdata.T',optMat,sprintf('example %d', matchRound));
%         fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
%             matchRound, optError,256*optError,fullError,256*fullError,100*overlapError);
%         fullErrors = [fullErrors, fullError];
    end
    matchRound = matchRound + 1;
end
nOtherMatch = matchRound - nFeatMatch - 1;

%% 7. Output the result.
% 似乎不用做什么。
% 不对，得把大小还原回去……
if ~useBigImages
    affines = RecoverResizedAffines(affines, matchConfigs, IResize, TResize);
end

fprintf('\n!!!!!!!!ENDING!!!!!!!!!!\n');
fprintf('nFeatMatch = %d, nOtherMatch = %d\n', nFeatMatch, nOtherMatch);
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
% 把用小图片找到的匹配还原成大图片的匹配

configs(:, 1:2) = configs(:, 1:2) / imgresize; % tx, ty
configs(:, 4:5) = configs(:, 4:5) / (imgresize/tplresize); % sx, sy

for i = 1:length(affs)
%     aff = affs(i);
    
    affmat = CreateAffineTransformation(configs(i,:));
    affs(i) = maketform('affine', affmat');
end

end

function [] = DisplayRestrictedSR(restrictedSearchRange, img_in)

figure, imshow(img_in);
for j=1:size(restrictedSearchRange,1)
rectangle('Position', [restrictedSearchRange(j,5) + size(img_in,2)/2, ...
    restrictedSearchRange(j,7) + size(img_in,1)/2, restrictedSearchRange(j,6) - restrictedSearchRange(j,5), ...
    restrictedSearchRange(j,8) - restrictedSearchRange(j,7) ]);
end

end
