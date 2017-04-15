function [affines, scores] = K_MultiMatch(I, T, params)
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
[BBSParams, matchParams] = CheckParams(params);

affines = [];
scores = [];

%% 1. Downsample _I_ and _T_ w.r.t. some rules
[I, T, resizeFactor] = AdjustSize(I, T, pz);
matchParams.resizeFactor = resizeFactor;

%% 2. Use BBS to approximately locate possible matches. Generate a mask _M_
[M] = BBSPreprocess(I, T, BBSParams);
szI = size(I);
szT = size(T);
%% 3. Generate configs for searching.
[initConfigs] = InitializeConfigs(M, matchParams, szI, szT);

%% 4. Search for the 1st match _A1_.
matchRound = 1;
fprintf('First match:\n');
[A, score, accept] = FindFirstMatch(I, T, M, initConfigs, matchParams);
% A is not a matrix, it is a struct!
if accept == 0
    fprintf('First match score (%.4f) too low, halt!\n', score);
    return;
else
    affines = [affines, A];
    scores = [scores, score];
end

%% 5. Set _M(A1(T))=0_ and search for the 2nd match, and so on.
curAccept = 1;
matchRound = matchRound + 1;
while(curAccept)
    fprintf('Match #%d:\n', matchRound);
    [M, initConfigs] = MaskNewMatch(M, initConfigs, A);
    [A, score, accept] = FindNextMatch(I, T, M, initConfigs, matchParams, scores); % scores can be used as threshold
    curAccept = accept;
    if(~curAccept)
        fprintf('Stop at match #%d!\n', matchRound);
    else
        affines = [affines, A];
        scores = [scores, score];
    end
end

%% 6. Output the result.
% 似乎不用做什么。

end

function [bbsparams, matchparams] = CheckParams(params)
% BBSParams: gamma, gz, rect(?)
% matchParams:
% epsilon, delta, photometricInvariance;
% searchRange: minTx, maxTx, minTy, maxTy, minRotation, maxRotation, minScale, maxScale

bbsparams = params.BBSParams;
matchparams = params.matchParams;

if(~isfield(bbsparams, 'gamma')); bbsparams.gamma = 2; end
if(~isfield(bbsparams, 'pz'));    bbsparams.pz = 3; end

if(~isfield(matchparams, 'epsilon')); matchparams.epsilon = 0.15; end
if(~isfield(matchparams, 'delta')); matchparams.delta = 0.25; end
if(~isfield(matchparams, 'photometricInvariance')); matchparams.photometricInvariance = 0; end
% searchRange will be initialized in InitializeConfigs()
if(~isfield(matchparams, 'searchRange')); matchparams.searchRange = []; end


end
