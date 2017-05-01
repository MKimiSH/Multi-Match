function [affines, scores, matchConfigs] = K_MultiMatchDBSCAN(img, tpl, params)
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


%% 1. Downsample _I_ and _T_ w.r.t. some rules
[img, tpl, IResize, TResize] = AdjustSize(img, tpl, BBSParams.pz);
matchParams.resizeFactor = IResize/TResize;

%% 2. Use BBS to approximately locate possible matches. Generate a mask _M_
[M] = BBSPreprocess(img, tpl, BBSParams);
img = MakeOdd(img);
tpl = MakeOdd(tpl);
M = MakeOdd(M);
szI = size(img);
szT = size(tpl);
%% 3. Generate configs for searching.
% M = ones(size(M)); % debug
[initConfigs, gridSize, matchParams] = InitializeConfigs(M, matchParams, szI, szT);

[goodConfigs, configClassList, numClasses] = FindGoodInitConfigs(img, tpl, M, initConfigs, matchParams);
%% 4. Search for the 1st match _A1_.
%% 5. Set _M(A1(T))=0_ and search for the 2nd match, and so on.
curAccept = 1;
matchRound = 1;
affines = [];
scores = zeros(numClasses, 1);
matchConfigs = zeros(numClasses, 6);
% matchRound = matchRound + 1;
A = [];
deltaFact = 1.511;
while(curAccept && matchRound<=5)
    fprintf('Match #%d, %d goodConfigs in this class: \n', matchRound, nnz(configClassList == matchRound));
    initGoodConfigs = goodConfigs(configClassList == matchRound, :);
    curConfigs = ExpandConfigsRandom(initGoodConfigs, matchParams.steps, 1, 80, deltaFact);
    curConfigs = BoundConfigsTrSc(curConfigs, matchParams.bounds);
    [M, initConfigs] = MaskNewMatch(M, tpl, initConfigs, A);
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
