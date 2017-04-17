function [mask] = BBSPreprocess(img, tpl, params)
% Use Best Buddy Similarity to approximately find the target search region
% I don't know to what extent can it handle rotation. Just a try.

[ih, iw, ~] = size(img);
[th, th, ~] = size(tpl);
szI = size(img);
szT = size(tpl);
gamma = params.gamma;
pz = params.pz;

%% calculate BBS
tic;
BBS = computeBBS(img, tpl, gamma, pz);
% interpolate likelihood map back to image size 
BBS = BBinterp(BBS, szT(1:2), pz, NaN); % wonder whether I should interpolate.
t = toc;
fprintf('BBS computed in %.2f sec (|I| = %dx%d , |T| = %dx%d)\n',t,szI(1:2),szT(1:2));

%% thresholding BBS value for a mask
meanBBS = mean(BBS(:), 'omitnan');
stdBBS = std(BBS(:), 'omitnan');
maxBBS = max(BBS(:), [], 'omitnan');
threshold = meanBBS + 1.5*stdBBS;
if(maxBBS < threshold)
    threshold = threshold - 0.5*stdBBS;
end

mask = BBS > threshold;
str = strel('disk',2); % 我也不知道要不要dilate
mask = imdilate(mask, str);
end