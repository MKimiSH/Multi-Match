dbstop if error

img = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\benchmark\targets\0005.jpg');
img = imresize(img, 2);
img = img(1:end-1, 1:end-1, :);

tpl = img(180:280, 270:370,:);

r1list = linspace(0,0,50);
r2list = linspace(0,0,50);
[R1, R2] = meshgrid(r1list, r2list);
rconfigs = [R1(:) R2(:)];

configs = zeros(length(rconfigs), 6);
configs(:,1) = 0;
configs(:,2) = 0;
configs(:,3) = rconfigs(:,2);
configs(:,4) = 1;
configs(:,5) = 1;
configs(:,6) = rconfigs(:,1);

[h1, w1, d1] = size(tpl);
[h2, w2, d2] = size(img);
r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);

[matrixConfigs_mex, insiders] = ...
    Configs2Affine_mex(configs',int32(h1), int32(w1), int32(h2), int32(w2), int32(r1x), int32(r1y), int32(r2x), int32(r2y));

inBoundaryInds = find(insiders);
matrixConfigs_mex = matrixConfigs_mex(:,inBoundaryInds);
origNumConfigs = size(configs,1);

configs = configs(inBoundaryInds,:);

mask = ones(h1, w1);
locs = find(mask);
numPoints = 1000;
ind = randi(length(locs), [1,numPoints]);
[ys,xs] = ind2sub(size(mask),locs(ind));
ys = ys';
xs = xs';

% ÎªÊ²Ã´»á¹Òµô°¡¡£¡£¡£
distances = EvaluateConfigsVectorized_mex(permute(tpl,[3,2,1]),permute(img,[3,2,1]),matrixConfigs_mex,int32(xs),int32(ys),int32(0));
