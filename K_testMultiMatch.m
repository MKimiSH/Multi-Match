function [] = K_testMultiMatch()

dbstop if error

AddPaths

img = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\3.jpg');
% img = rgb2gray(img);
img = im2double(img);
% img = imrotate(img, 90); % poor result with rotation.
tpl = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\t3.png');
% tpl = rgb2gray(tpl);
tpl = im2double(tpl);

% img = rgb2gray(img); tpl = rgb2gray(tpl);

% MultiMatch run
params = [];
tic;
% [affines, scores, mConfigs] = K_MultiMatchRound(img, tpl, params);
[affines, scores, mConfigs] = K_MultiMatchDBSCAN(img, tpl, params);
tMatch = toc;

% bestTransMat = affines(1).tdata.T';
optMat = [];
% Visualize results, one affine per image
for i=1:length(affines)
    bestTransMat = affines(i).tdata.T';
    [optError,fullError,overlapError] = MatchingResult(tpl,img,bestTransMat,optMat,sprintf('example %d', i));
    fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        i, optError,256*optError,fullError,256*fullError,100*overlapError);
    fprintf('example %d: finished\n\n', i);
end

fprintf('Match Time (%d matches): %.4f seconds.\n', length(affines), tMatch);

keyboard;
end