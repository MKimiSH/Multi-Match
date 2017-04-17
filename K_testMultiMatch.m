dbstop if error

AddPaths

img = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\1.jpg');
% img = rgb2gray(img);
img = im2double(img);
tpl = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\t1.png');
% tpl = rgb2gray(tpl);
tpl = im2double(tpl);

img = rgb2gray(img); tpl = rgb2gray(tpl);

% MultiMatch run
params = [];
[affines, scores, mConfigs] = K_MultiMatch(img, tpl, params);

bestTransMat = affines(1).tdata.T';
optMat = [];
% Visualize result
[optError,fullError,overlapError] = MatchingResult(tpl,img,bestTransMat,optMat,'example 1');
fprintf('example 1 - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
    optError,256*optError,fullError,256*fullError,100*overlapError);
fprintf('example 1: finished\n\n');

bestTransMat = affines(2).tdata.T';
[optError,fullError,overlapError] = MatchingResult(tpl,img,bestTransMat,optMat,'example 2');
fprintf('example 2 - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
    optError,256*optError,fullError,256*fullError,100*overlapError);
fprintf('example 2: finished\n\n');