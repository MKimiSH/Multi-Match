dbstop if error
AddPaths

img = imread('testimg\original\class1\01.jpg');
% img = rgb2gray(img);
img = im2double(img);
img = MakeOdd(img);
tpl = imread('testimg\templates\class1\5.png');
% tpl = rgb2gray(tpl);
tpl = im2double(tpl);
tpl = MakeOdd(tpl);

% FastMatch run
[bestConfig,bestTransMat,sampledError] = FastMatch(tpl,img);

optMat = [];
% Visualize result
[optError,fullError,overlapError] = MatchingResult(tpl,img,bestTransMat,optMat,'example 1');

fprintf('example 1 - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
    optError,256*optError,fullError,256*fullError,100*overlapError);
fprintf('example 1: finished\n\n');