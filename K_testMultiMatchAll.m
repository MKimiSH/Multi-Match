function [] = K_testMultiMatchAll()
% ALLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL!

dbstop if error
clc;
close all;

AddPaths

testpairs = [1,1; 1,2; 2,2; 3,3; 4 4;5 5; 6 6; 7 7; 7 10; 8 6; 
            9 6; 10 6; 12 8; 12 11; 13 9; 14 7; 15 8; 15 11;
            16 12; 17 13; 18 14; 19 15; 20 16; 21 17; 22 18; 23 19; 24 20; 25 21];

img = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\24.jpg');
% img = rgb2gray(img);
img = im2double(img);
% img = imrotate(img, 90); % poor result with rotation.
tpl = imread('C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\t20.png');
% tpl = rgb2gray(tpl);
tpl = im2double(tpl);

% img = rgb2gray(img); tpl = rgb2gray(tpl);

% MultiMatch run
params = [];
tic;
% [affines, scores, mConfigs] = K_MultiMatchRound(img, tpl, params);
[affines, scores, mConfigs, timeList, fullErrors] = K_MultiMatchFeaturePoints(img, tpl, params);
tMatch = toc;

% bestTransMat = affines(1).tdata.T';
optMat = [];
% Visualize results, one affine per image
% for i=1:length(affines)
%     bestTransMat = affines(i).tdata.T';
%     [optError,fullError,overlapError] = MatchingResult(tpl,img,bestTransMat,optMat,sprintf('example %d', i));
%     fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
%         i, optError,256*optError,fullError,256*fullError,100*overlapError);
%     fprintf('example %d: finished\n\n', i);
% end

fprintf('Match Time (%d matches): %.4f seconds.\n', length(affines), tMatch);

keyboard;
close all;
end