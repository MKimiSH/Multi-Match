function [configs, gridSize, params] = InitializeConfigs(mask, params, szI, szT)
% Generate initial configs for searching
% _configs_: tx, ty, r2, sx, sy, r1 (like FastMatch)
% 傻了，人家是以图像中心(r2x, r2y)为原点的
% tpl也是以tpl中心为原点的

%% init searchRange
searchRange = params.searchRange;
resizeFactor = params.resizeFactor;
r1x = 0.5*(szT(2)-1);
r1y = 0.5*(szT(1)-1);
r2x = 0.5*(szI(2)-1);
r2y = 0.5*(szI(1)-1);
if ~isfield(searchRange,'minScale'), searchRange.minScale = 0.5 * resizeFactor; end
if ~isfield(searchRange,'maxScale'), searchRange.maxScale = 2 * resizeFactor; end
if ~isfield(searchRange,'minRotation'), searchRange.minRotation = -pi; end
if ~isfield(searchRange,'maxRotation'), searchRange.maxRotation = pi; end
if ~isfield(searchRange,'minTx'), searchRange.minTx = -(r2x-r1x*searchRange.minScale); end
if ~isfield(searchRange,'maxTx'), searchRange.maxTx = r2x-r1x*searchRange.minScale; end
if ~isfield(searchRange,'minTy'), searchRange.minTy = -(r2y-r1y*searchRange.minScale); end
if ~isfield(searchRange,'maxTy'), searchRange.maxTy = r2y-r1y*searchRange.minScale; end

% copy params
delta = params.delta;
minScale = searchRange.minScale;
maxScale = searchRange.maxScale;
minRotation = searchRange.minRotation;
maxRotation = searchRange.maxRotation;
minTx = max(searchRange.minTx,-(r2x-r1x*minScale));
maxTx = min(searchRange.maxTx,r2x-r1x*minScale);
minTy = max(searchRange.minTy,-(r2y-r1y*minScale));
maxTy = min(searchRange.maxTy,r2y-r1y*minScale);

%% parametrize the initial grid
[bounds,steps] = GenerateGrid(szT(2),szT(1),delta,minTx,maxTx,minTy,maxTy,minRotation,maxRotation,minScale,maxScale);
inlierRatio = nnz(mask)/numel(mask);
% steps.tx = steps.tx * sqrt(inlierRatio);
% steps.ty = steps.ty * sqrt(inlierRatio);

params.searchRange = searchRange;
params.bounds = bounds;
params.steps = steps;

%% use _mask_ to restrict _tx_ and _ty_

tx_steps = bounds.tx(1) : steps.tx : (bounds.tx(2) + 0.5*steps.tx); % 'pad' at end of range with an extra sample
ty_steps = bounds.ty(1) : steps.ty : (bounds.ty(2) + 0.5*steps.ty); % 'pad' at end of range with an extra sample
ntx_steps = length(tx_steps);
nty_steps = length(ty_steps);

txtymask = zeros(nty_steps, ntx_steps);
% 先写个for，上来就优化太烦了。。
% [txgrid, tygrid] = meshgrid(tx_steps, ty_steps);
% txgrid = round(txgrid);
% tygrid = round(tygrid);

mv_tx_steps = tx_steps + r2x; % move to center
mv_ty_steps = ty_steps + r2y;

tyWithinArea = mv_ty_steps>0.5 & mv_ty_steps<=szI(1);
txWithinArea = mv_tx_steps>0.5 & mv_tx_steps<=szI(2);
for i = 1:nty_steps
    if tyWithinArea(i)
        for j = 1:ntx_steps
            if txWithinArea(j) && mask(round(mv_ty_steps(i)), round(mv_tx_steps(j)))
                txtymask(i,j) = 1;
            end
        end
    end
end
txty_steps = zeros(nnz(txtymask), 2);
[row, col] = find(txtymask==1);
txty_steps(:,1) = tx_steps(col);
txty_steps(:,2) = ty_steps(row);
ntxty_steps = nnz(txtymask);

%% Prepare other configs, r2, sx, sy, r1
% Rotations ignore the user selected range here - it is handles in FindBestTransformation
r_steps = -pi : steps.r : pi; % no padding since it is a cyclic range
s_steps = bounds.s(1) : steps.s : bounds.s(2) + 0.5*steps.s; % 'pad' at end of range with an extra sample
if (steps.s == 0)
    s_steps = bounds.s(1);
end
ns_steps = length(s_steps);
nr_steps = length(r_steps);

% second rotation is a special case (can be limited to a single quartile)
quartile1_r_steps = r_steps(r_steps < -pi/2 + steps.r/2);
NR2_steps = length(quartile1_r_steps);

% gridsize
gridSize = ntxty_steps*(ns_steps^2)*(nr_steps*NR2_steps);

%% Create _configs_
configs = zeros(gridSize, 6);

% MATLAB 好强啊！
[T, R2, SX, SY, R1] = ndgrid(1:ntxty_steps, quartile1_r_steps, s_steps, s_steps, r_steps);
configs(:,2:6) = [T(:), R2(:), SX(:), SY(:), R1(:)];
configs(:,1:2) = txty_steps(configs(:,2), :); % 这尴尬的indexing...

% % 毅种循环
% txty_interval = gridSize / ntxty_steps;
% r2_interval = txty_interval / NR2_steps;
% sx_interval = r2_interval / ns_steps;
% sy_interval = sx_interval / ns_steps;
% r1_interval = 1; % = sy_interval / nr_steps
% % 跳跃步长
% txty_jump = 1;
% r2_jump = txty_jump * ntxty_steps;
% sx_jump = r2_jump * NR2_steps;
% sy_jump = sx_jump * ns_steps;
% r1_jump = sy_jump * ns_steps;
% 
% for txty_ind = 1:ntxty_steps
%     configs((txty_ind-1)*txty_interval + 1:txty_jump:txty_ind*txty_interval, 1) = txty_steps(txty_ind, 1);
%     configs((txty_ind-1)*txty_interval + 1:txty_jump:txty_ind*txty_interval, 2) = txty_steps(txty_ind, 2);
% end
% for r2_ind = 1:NR2_steps
%     configs((r2_ind-1)*r2_interval + 1:r2_jump:r2_ind*r2_interval, 3) = r2_steps(r2_ind);
% end
% for sx_ind = 1:ns_steps
%     configs((sx_ind-1)*sx_interval + 1:sx_jump:sx_ind*sx_interval, 4) = sx_steps(sx_ind);
% end
% for sy_ind = 1:ny_steps
%     configs((sy_ind-1)*sy_interval + 1:sy_jump:sy_ind*sy_interval, 5) = sy_steps(sy_ind);
% end
% for r1_ind = 1:nr_steps
%     configs((r1_ind-1)*r1_interval + 1:r1_jump:r1_ind*r1_interval, 6) = r1_steps(r1_ind);
% end 
end