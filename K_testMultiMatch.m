function [] = K_testMultiMatch()

dbstop if error
clc;
close all;

AddPaths

testpairs = [1,1; 1,2; 2,2; 3,3; 4 4;5 5; 6 6; 7 7; 7 10; 8 6; 
            9 6; 10 6; 12 8; 12 11; 13 9; 14 22; 15 8; 15 11;
            16 12; 17 13; 18 14; 19 15; 20 16; 21 17; 22 18; 23 19; 24 20; 25 21];
matchNums = [10; 4; 20; 10; 6; 5; 7; 14; 2; 4; 6; 10; 13; 5; 15; 8; 3; 25; 7;
            2; 9; 7; 11; 11; 4; 4; 3; 9];   

        
% set(0,'DefaultFigureVisible', 'off')
for useBBS = [1 0]
    for useSURF = [1 0]
        for testNum = 1:length(testpairs)
            clc;
            close all;
            
            % MultiMatch run
            params = [];
            params.nMatches = matchNums(testNum);
            params.useBBS = useBBS;
            params.useSURF = useSURF;
            
%             imgNum = testpairs(testNum, 1);
%             tplNum = testpairs(testNum, 2);
            imgNum = 15; tplNum = 11;
            dirName = ['results4\', num2str(imgNum), '_t', num2str(tplNum), '_BBS', num2str(params.useBBS), '_SURF', num2str(params.useSURF), '\']
            if(exist([dirName, 'res.mat'], 'file'))
                continue;
            end
            if(imgNum==6), continue; end
            
            mkdir(dirName);
            img = imread(['C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\', num2str(imgNum), '.jpg']);
            if(size(img,1)>3000 || size(img,2)>3000), continue; end
            % img = rgb2gray(img);
            img = im2double(img);
            % img = imrotate(img, 90); % poor result with rotation.
            tpl = imread(['C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\t',num2str(tplNum), '.png']);
            % tpl = rgb2gray(tpl);
            tpl = im2double(tpl);
            
            
            % [affines, scores, mConfigs] = K_MultiMatchRound(img, tpl, params);
            [affines, scores, mConfigs, timeList, fullErrors, nFeatureMatch, nOtherMatch] = K_MultiMatchFeaturePoints(img, tpl, params);
            timeList
            fmt = sum(timeList.featureMatchTime)
            omt = sum(timeList.otherMatchTime)
            % figure, plot(scores)
            optMat = [];
            % Visualize results, one affine per image
            for i=1:length(affines)
                bestTransMat = affines(i).tdata.T';
                [optError,fullError,overlapError, h] = MatchingResult(tpl,img,bestTransMat,optMat,sprintf('example %d', i));
                fprintf('example %d - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
                    i, optError,256*optError,fullError,256*fullError,100*overlapError);
                fullErrors = [fullErrors, fullError];
                fprintf('example %d: finished\n\n', i);
                saveas(h, [dirName, sprintf('match%02d.png', i)]);
            end
            % fprintf('Match Time (%d matches): %.4f seconds.\n', length(affines), tMatch);
            save([dirName, 'res.mat'], 'affines', 'scores', 'mConfigs', 'timeList', 'fullErrors', 'nFeatureMatch', 'nOtherMatch');
        end
    end
end
set(0,'DefaultFigureVisible', 'on')
% keyboard;
close all;
end