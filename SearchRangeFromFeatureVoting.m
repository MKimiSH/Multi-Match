function [restrictedSearchRanges, resM, scaleRotRange, groups, srFromVote, time] = SearchRangeFromFeatureVoting(img, tpl, M)
% ��img��tpl���ҵ�ƥ���SURF����
% ������λ����Ҫ��M��
% Ȼ��ƥ��������������ͶƱ
% ͶƱת��Ϊ������Χ��
% ������Χ������һ���ṹ�壬����һ���ṹ�塣

tic;
restrictedSearchRanges = [];

votes = [];

[ih, iw, ~] = size(img);
img_gray = rgb2gray(img);
tpl_gray = rgb2gray(tpl);

%% Detect and match SURF features in both images
% use M to select inliers;
resM = imresize(M, [ih, iw]);
% blur img_gray!
dilatedM = imresize(imdilate(M, strel('square',5)), [ih, iw]);
% img_blurredgray = imfilter(img_gray, fspecial('gaussian', [30 30], 3));
% img_combgray = img_blurredgray;
% img_combgray(dilatedM) = img_gray(dilatedM);

points_img = detectSURFFeatures(img_gray, 'MetricThreshold', 800);
points_tpl = detectSURFFeatures(tpl_gray, 'MetricThreshold', 800);

% �ҵ���Mask��Χ�ڵ�������
isInBound = (round(points_img.Location(:,1)) > 0 & round(points_img.Location(:,1)) <= iw ...
    & round(points_img.Location(:,2)) > 0 & round(points_img.Location(:,2)) <= ih);
locations = points_img.Location(isInBound,:);
roundedLocations = round(locations);
indLocations = sub2ind(size(resM), roundedLocations(:,2), roundedLocations(:,1));
isInM = dilatedM(indLocations) > 0;
inBoundIdx = find(isInBound);
inMIdx = inBoundIdx(isInM);

goodPoints_img = points_img(inMIdx);
% goodPoints_img = points_img;
goodPoints_tpl = points_tpl;

% extract features
[feat_img, vpoint_img] = extractFeatures(img_gray, goodPoints_img);
[feat_tpl, vpoint_tpl] = extractFeatures(tpl_gray, goodPoints_tpl);

indexPairs = matchFeatures(feat_img, feat_tpl);

if isempty(indexPairs)
    scaleRotRange = [];
    groups = [];
    srFromVote = [];
    time = toc;
    fprintf('No matching features found.\n');
%     keyboard;
    return;
end
% % centering feat_tpl
% vpoint_tpl_ctrd = vpoint_tpl;
% [th, tw, ~] = size(tpl);
% for i=1:length(vpoint_tpl)
%     vpoint_tpl_ctrd(i).Location = vpoint_tpl_ctrd(i).Location - [(tw-1)/2, (th-1)/2];
% end

%% Calculate vote by matched features
[th, tw, ~] = size(tpl);
nPairs = size(indexPairs, 1);
votes = zeros(nPairs, 4); % tx, ty, s, r
for i=1:nPairs
    fim = vpoint_img(indexPairs(i,1));
    ftp = vpoint_tpl(indexPairs(i,2));
    vote_r = fim.Orientation - ftp.Orientation; % rotation
    vote_s = fim.Scale / ftp.Scale; % scale
    rotMat = [cos(vote_r) -sin(vote_r); sin(vote_r) cos(vote_r)];
    vote_t = fim.Location' - vote_s * (rotMat*(ftp.Location' - [(tw-1)/2, (th-1)/2]')); % Location need to be a column vector
    vote_r = mod(vote_r + 11*pi, 2*pi) - pi; % vote_r \in [-pi, pi)
    votes(i,:) = [vote_t(1), vote_t(2), vote_s, vote_r];
end
%% Calc. search ranges from votes.
srFromVote = zeros(nPairs, 8); % mins, maxs, minr, maxr, mintx, maxtx, minty, maxty
scaleFact = 0.2;
translateFact = 15;
rotationFact = 0.25;

% vote: tx, ty, s, r
mins = votes(:,3) - scaleFact;
maxs = votes(:,3) + scaleFact;
minr = votes(:,4) - rotationFact;
maxr = votes(:,4) + rotationFact;
mintx = votes(:,1) - votes(:,3)*translateFact;
maxtx = votes(:,1) + votes(:,3)*translateFact;
minty = votes(:,2) - votes(:,3)*translateFact;
maxty = votes(:,2) + votes(:,3)*translateFact; % ��д�������ΰɣ�
srFromVote = [mins, maxs, minr, maxr, mintx, maxtx, minty, maxty];

%% Group Search Ranges
nbMat = zeros(nPairs, nPairs); % ���ο��า�Ǿ���Ϊ������ 
% ��ֻ��ƽ�ưɣ�
for i=1:nPairs
    cur = srFromVote(i,:);
%     rect = cur(5:8);
    nbMat(i,:) = ~(srFromVote(:,5) > cur(6) | srFromVote(:,6) < cur(5) ...
            | srFromVote(:,7) > cur(8) | srFromVote(:,8) <cur(7));
end

% ���¾��࣬��������ȶ������Ҫ�ȶ�Ҫ��DBSCAN��
groups = zeros(nPairs,1);
nGroups = 0;
curGroup = 1;
for i=1:nPairs
    nbs = nbMat(i,:);
    nbs(1:i) = 0;
    if groups(i)>0
        groups(nbs>0) = groups(i);
    else
        if any(groups(nbs>0))
            avalGroups = groups>0;
            avalGroups = avalGroups(nbs>0);
            avalGroups = groups(avalGroups);
            groups(i) = avalGroups(1);
        else
            groups(i) = curGroup;
            curGroup = curGroup + 1;
        end
        groups(nbs>0) = groups(i);
    end
end
nGroups = curGroup - 1;
fprintf('......nGroups = %d\n', nGroups);

% ��ÿ�������������Χ
restrictedSearchRanges = zeros(nGroups, 8);
for g = 1:nGroups
    cands = srFromVote(groups==g,:); %candidates
    mins = min(cands(:,1));
    maxs = max(cands(:,2));
    minr = min(cands(:,3));
    maxr = max(cands(:,4));
    mintx = min(cands(:,5)) - (iw-1)/2;
    maxtx = max(cands(:,6)) - (iw-1)/2;
    minty = min(cands(:,7)) - (ih-1)/2;
    maxty = max(cands(:,8)) - (ih-1)/2;
    restrictedSearchRanges(g,:) = [mins, maxs, minr, maxr, mintx, maxtx, minty, maxty];
end

% �����scale��rotation������
allmins = min(restrictedSearchRanges(:,1));
allmaxs = max(restrictedSearchRanges(:,2));
allminr = min(restrictedSearchRanges(:,3));
allmaxr = max(restrictedSearchRanges(:,4));
scaleRotRange = [allmins, allmaxs, allminr, allmaxr];

time = toc;
end