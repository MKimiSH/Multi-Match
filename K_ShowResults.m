function [] = K_ShowResults()

testpairs = [1,1; 1,2; 2,2; 3,3; 4 4;5 5; 6 6; 7 7; 7 10; 8 6; 
            9 6; 10 6; 12 8; 12 11; 13 9; 14 22; 15 8; 15 11;
            16 12; 17 13; 18 14; 19 15; 20 16; 21 17; 22 18; 23 19; 24 20; 25 21];
matchNums = [10; 4; 20; 10; 6; 5; 7; 14; 2; 4; 6; 10; 13; 5; 15; 8; 3; 25; 7;
            2; 9; 7; 11; 11; 4; 4; 3; 9];   
        
showImage = 1;
for useBBS = [1 0]
    for useSURF = [1 0]
        for testNum = 1%:length(testpairs)
%             imgNum = testpairs(testNum, 1);
%             tplNum = testpairs(testNum, 2);
            imgNum = 18; tplNum = 14;
            dirName = ['results3\', num2str(imgNum), '_t', num2str(tplNum), '_BBS', num2str(useBBS), '_SURF', num2str(useSURF), '\'];
            if(~exist([dirName, 'res.mat'], 'file'))
                fprintf('0\t0\t0\t0\n');
                continue;
            end
%             if(imgNum==6); continue; end
            load([dirName, 'res.mat']);
            fmt = sum(timeList.featureMatchTime);
            omt = sum(timeList.otherMatchTime);
%             fprintf('%d.jpg -- t%d.png\n', imgNum, tplNum);
%             fprintf('useBBS=%d, useSURF=%d\n', useBBS, useSURF);
            if exist('nFeatureMatch', 'var')
%                 fprintf('nFeatureMatch = %d, nOtherMatch = %d\n', nFeatureMatch, nOtherMatch);
            end
            if(useBBS)
                fprintf('%.4f\t', timeList.BBSTime);
            end
            if(useSURF)
                fprintf('%.4f\t%.4f\t', timeList.featureTime, fmt);
            end
            fprintf('%.4f\n', omt);
            %             fprintf('-----------------------------------------------\n');
            if(showImage)
                gt = load(['C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\matches\',  num2str(imgNum), '_t', num2str(tplNum), '.mat']);
                gtConfigs = gt.configs;
                img = imread(['C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\', num2str(imgNum), '.jpg']);
                img = im2double(img);
                tpl = imread(['C:\M.K.S.H\Study and study\Scientific work\Affine matching\all_imgs\dataset_K\t',num2str(tplNum), '.png']);
                tpl = im2double(tpl);
                MultipleMatchesInOneImage(img, tpl, mConfigs, gtConfigs);
            end
        end
    end
end
end