function [configs] = InitializeConfigs(mask, params, szI, szT)
% Generate initial configs for searching
% _configs_: tx, ty, r2, sx, sy, r1 (like FastMatch)

%% init searchRange
searchRange = params.searchRange;
r1x = 0.5*(szI(2)-1);
r1y = 0.5*(szI(1)-1);
r2x = 0.5*(szT(2)-1);
r2y = 0.5*(szT(1)-1);
if ~isfield(searchRange,'minScale'), searchRange.minScale = 0.5; end
if ~isfield(searchRange,'maxScale'), searchRange.maxScale = 2; end
if ~isfield(searchRange,'minRotation'), searchRange.minRotation = -pi; end
if ~isfield(searchRange,'maxRotation'), searchRange.maxRotation = pi; end
if ~isfield(searchRange,'minTx'), searchRange.minTx = -(r2x-r1x*searchRange.minScale); end
if ~isfield(searchRange,'maxTx'), searchRange.maxTx = r2x-r1x*searchRange.minScale; end
if ~isfield(searchRange,'minTy'), searchRange.minTy = -(r2y-r1y*searchRange.minScale); end
if ~isfield(searchRange,'maxTy'), searchRange.maxTy = r2y-r1y*searchRange.minScale; end

% copy params
minScale = searchRange.minScale;
maxScale = searchRange.maxScale;
minRotation = searchRange.minRotation;
maxRotation = searchRange.maxRotation;
minTx = max(searchRange.minTx,-(r2x-r1x*minScale));
maxTx = min(searchRange.maxTx,r2x-r1x*minScale);
minTy = max(searchRange.minTy,-(r2y-r1y*minScale));
maxTy = min(searchRange.maxTy,r2y-r1y*minScale);

%% parametrize the initial grid
[bounds,steps] = GenerateGrid(w1,h1,delta,minTx,maxTx,minTy,maxTy,minRotation,maxRotation,minScale,maxScale);
inlierRatio = nnz(mask)/numel(mask);
steps.tx = steps.tx * sqrt(inlierRatio);
steps.ty = steps.ty * sqrt(inlierRatio);

%% use _mask_ to restrict _tx_ and _ty_
inliers = find(mask==1);


end