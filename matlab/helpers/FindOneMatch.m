function [aff, score, config, accept] = FindOneMatch(img, tpl, mask, configs, params, scores)
% Use FastMatch-like method to find the first match!
%% verify input image types
I2 = img;
I1 = tpl;
szI = size(img);
szT = size(tpl);
if ( ~strcmp(class(I1),'double') || ~strcmp(class(I2),'double')) %#ok<STISA>
    error('FastMatch: I1 and I2 should both be of class ''double'' (in the range [0,1])');
end

templateMask = ones(size(I1, 1), size(I1, 2));
if ((size(templateMask,1) ~= size(I1,1)) || (size(templateMask,2) ~= size(I1,2)) )
    error('FastMatch: Template mask not same size as template');
end

isGrayscale = (size(I1,3)==1);

if ~isGrayscale
%     error('FastMatch: Only grayscale images are currently supported');
end
%% take out params
if ~exist('scores', 'var')
    scores = [];
    stage = 'first';
else
    stage = 'next';
end
epsilon = params.epsilon;
delta = params.delta;
photometricInvariance = params.photometricInvariance;
bounds = params.bounds;
steps = params.steps;

%% blur in main loop - this reduces the total-variation and gives better results
if (isGrayscale)
    origI1 = I1;
    origI2 = I2;
end
[h1,w1,d] = size(I1);

%% generate Theta(1/eps^2) random points (and fresh ones each iteration later on)
numPoints = round(10/epsilon^2);
[xs, ys] = getPixelSamplePatch(templateMask, numPoints); %改

%% generate the Net
% [configs,gridSize] = CreateListOfConfigs(bounds,steps); %已经有了

if (size(configs,1) > 71000000)
        error('more than 35 million configs!');
end


%% main loop

deltaFact = 1.511;
level = 0;
bestDists = [];
perRoundNumConfigs = [];
perRoundNumGoodConfigs = [];
perRoundOrig_percentage = [];
bestGridVec = [];
newDelta = delta;
totTime = 0;
while (1)
        level = level + 1;
        
        if (isGrayscale) % slightly blur to reduce total-variation
            blur_sigma = 1.5+0.5/deltaFact^(level-1); % 2;
            blur_size = ceil(4 * blur_sigma);
            params.blur_kernel  = fspecial('gaussian', blur_size, blur_sigma);
            
            I1 = imfilter(origI1,params.blur_kernel,'symmetric');
            I2 = imfilter(origI2,params.blur_kernel,'symmetric');
        end
        
        [h2,w2,d2] = size(I2);
        
        r1x = 0.5*(w1-1);
        r1y = 0.5*(h1-1);
        r2x = 0.5*(w2-1);
        r2y = 0.5*(h2-1);
        
        % 0] if limited rotation range - filter out illegal rotations
        if (bounds.r(1)>-pi || bounds.r(2)<pi)
            minRot = bounds.r(1);
            maxRot = bounds.r(2);
            % total rotation in the range [0,2*pi]
            totalRots = mod(configs(:,3)+configs(:,6),2*pi);
            % total rotation in the range [-pi,pi]
            totalRots(totalRots>pi) = totalRots(totalRots>pi) - 2*pi;
            % filtering
            configs = configs(totalRots>=minRot & totalRots<=maxRot,:);
        end
        
        % 1] translate config vectors to matrix form
        Configs2AffineMEX = tic;
        fprintf('----- Configs2Affine, with %d configs -----\n',size(configs,1));
%         configs = [0 0 0 1 1 0]; % [tx,ty,r2,sx,sy,r1]
        [matrixConfigs_mex, insiders] = ...
                Configs2Affine_mex(configs',int32(h1), int32(w1), int32(h2), int32(w2), int32(r1x), int32(r1y), int32(r2x), int32(r2y));
        
        inBoundaryInds = find(insiders);
        matrixConfigs_mex = matrixConfigs_mex(:,inBoundaryInds);
        origNumConfigs = size(configs,1);
        
        configs = configs(inBoundaryInds,:);
        Configs2Affine_mex_time = toc(Configs2AffineMEX);
        
        % 2] evaluate all configurations
        EvaluateConfigsMEX = tic;
        
        if (isGrayscale)
            distances = EvaluateConfigs_mex(I1',I2',matrixConfigs_mex,int32(xs),int32(ys),int32(photometricInvariance));
            fprintf('----- Evaluate Configs grayscale, with %d configs -----\n',size(configs,1));
        else
            distances = EvaluateConfigsVectorizedSAD_mex(permute(I1,[3,2,1]),permute(I2,[3,2,1]),matrixConfigs_mex,int32(xs),int32(ys),int32(photometricInvariance));
            fprintf('----- Evaluate Configs vectorized, with %d configs -----\n',size(configs,1));
        end
        
        EvaluateConfigs_mex_time = toc(EvaluateConfigsMEX);
        
        totTime = totTime + Configs2Affine_mex_time + EvaluateConfigs_mex_time;
        
        [bestDist,ind] = min(distances);
        bestConfig = configs(ind,:);
        bestTransMat = CreateAffineTransformation(configs(ind,:));        
        
        
        % 3] choose the 'surviving' configs and delta for next round

        [goodConfigs,tooHighPercentage,extremelyHighPercentage,veryLowPercentage,orig_percentage,thresh] = ...
            GetGoodConfigsByDistance(configs,bestDist,newDelta,distances,bestGridVec);
        
        numGoodConfigs = size(goodConfigs,1);

        fprintf('$$$ bestDist = %.3f\n',bestDist);
        fprintf('$$ numGoodConfigs: %d (out of %d), orig percentage: %.4f, bestDist: %.4f, thresh: %.4f\n',...
            size(goodConfigs,1), size(configs,1), orig_percentage, bestDist, thresh);
        
        % collect round stats
        bestDists(level) = bestDist; %#ok<AGROW>
        perRoundNumConfigs(level) = origNumConfigs; %#ok<AGROW>
        perRoundNumGoodConfigs(level) = numGoodConfigs; %#ok<AGROW>
        perRoundOrig_percentage(level) = orig_percentage; %#ok<AGROW>

        % 4] break conditions of Branch-and-Bound
        clear conditions
        conditions(1) = (bestDist < 0.005); % good enough 1
        conditions(2) = (level > 5) && (bestDist < 0.01); % good enough 2
        conditions(3) = (level >= 20); % enough levels
        conditions(4) = ((level > 3) && (bestDist > mean(bestDists(level-3:level-1))*0.97)); % no improvement in last 3 rounds
        conditions(5) = ((level > 2) && (numGoodConfigs>1000) && extremelyHighPercentage ); % too high expansion rate
        conditions(6) = ((level > 3) && (numGoodConfigs>1000) && (numGoodConfigs>50*min(perRoundNumGoodConfigs))); % a deterioration in the focus

        if any(conditions)
            fprintf('breaking BnB at level %d due to conditions: %s\n', level, num2str(find(conditions)))
            fprintf('best distances by round:  %s\n', num2str(bestDists,'  %.4f'))
            fprintf('num configs per round:    %s\n', num2str(round(perRoundNumConfigs/1000),'  %.5dK'))
            fprintf('# good configs per round: %s\n', num2str(round(perRoundNumGoodConfigs),'  %.6d'))
            fprintf('percentage to expand:     %s\n', num2str(perRoundOrig_percentage,'  %.4f'))
            break
        end
        
        
        % 6] debug: visualize on histogram
        
        
        % 7] expand 'surviving' configs for next round
             % ('restart' = [with smaller delta] if not too many configs and not too high percentage of configs to expand)
        if (~veryLowPercentage && ... % && ...
            ( (tooHighPercentage && (bestDist > 0.1) && ((level==1) && (origNumConfigs < 7.5*10^6)) ) || ...
              (                     (bestDist > 0.15)  && ((level==1) && (origNumConfigs <   5*10^6)) ) ) )
                fact = 0.9;
                fprintf('##### RESTARTING!!! changing from delta: %.3f, to delta: %.3f\n', newDelta, newDelta*fact);
                newDelta = newDelta*fact;
                level = 0;
                steps.tx = fact*steps.tx;
                steps.ty = fact*steps.ty;
                steps.r = fact*steps.r;
                steps.s = fact*steps.s;
                [configs,gridSize] = InitializeConfigs(mask, params, szI, szT); % 原来是CreateListOfConfigs()
        else
                prevDelta = newDelta;
                newDelta = newDelta/deltaFact;
                fprintf('##### CONTINUING!!! prevDelta = %.3f,  newDelta = %.3f \n',prevDelta,newDelta);
                
                % expand the good configs
                expandType = 'randomExpansion'; %  'fullExpansion'; %  'deltaGrid'; %
                switch expandType
                        case 'randomExpansion'
                                expandedConfigs = ExpandConfigsRandom(goodConfigs,steps,level,80,deltaFact);
                        case 'fullExpansion'
                                expandedConfigs = ExpandConfigsFull(goodConfigs,steps,level,deltaFact);
                end
                configs = [goodConfigs ; expandedConfigs];
                configs = BoundConfigsTrSc(configs, bounds);
                configs = MaskOutlierConfigs(configs, mask);
        end
        
        
        %     configs = unique(configs,'rows'); % REMOVED THIS - IT IS WORTHWHILE
        
        fprintf('***\n');
        fprintf('*** level %d:|goodConfigs| = %d, |expandedConfigs| = %d\n',level,numGoodConfigs,size(configs,1));
        fprintf('***\n');
        
        
        % 8] refresh random points
        [xs, ys] = getPixelSamplePatch(templateMask, numPoints);
end

%% debug error
if isGrayscale
    [xs,ys] = meshgrid(1:w1,1:h1);
    xs = xs(:)';
    ys = ys(:)';
    BD_full = EvaluateConfigs_mex(I1',I2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    BD_full_orig = EvaluateConfigs_mex(origI1',origI2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    fprintf('bestDist: %.4f, BD_full: %.4f, BD_full_orig: %.4f\n',bestDist,BD_full,BD_full_orig);
else
    [xs,ys] = meshgrid(1:w1,1:h1);
    xs = xs(:)';
    ys = ys(:)';
    BD_full = EvaluateConfigsVectorized_mex(permute(I1,[3,2,1]), permute(I2,[3,2,1]), bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    % BD_full_orig = EvaluateConfigsVectorized_mex(origI1',origI2',bestTransMat([1 4 7 2 5 8])',int32(xs),int32(ys),int32(photometricInvariance));
    fprintf('Color image: bestDist: %.4f, BD_full: %.4f\n',bestDist,BD_full);
end

%% for output
sampledError = bestDist;
score = sampledError;
aff = maketform('affine', bestTransMat');
config = bestConfig;
if strcmp(stage, 'first')
    accept = score < 0.2;
else
    accept = score < 0.4; % 先这样！
end

return
end



function [res,i] = IsMemberApprox(A,row,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;
for i = 1 : size(A,1)
        if (norm(A(i,:)-row) < err)
                res = 1;
                return
        end
end
end


function [goodConfigs,tooHighPercentage,extremelyHighPercentage,veryLowPercentage,orig_percentage,thresh] = ...
    GetGoodConfigsByDistance(configs,bestDist,newDelta,distances,bestGridVec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targetNum = 20000;
% thresh = bestDist + newDelta/3;
thresh = bestDist + GetThreshPerDelta(newDelta);
goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
numGoodConfigs = size(goodConfigs,1);
orig_percentage = numGoodConfigs/size(configs,1);

% too many good configs - reducing threshold
while (numGoodConfigs > 27000)
        thresh = thresh * 0.99;
        goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
        numGoodConfigs = size(goodConfigs,1);
end

if (isempty(goodConfigs))
         thresh = min(distances);
        goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
        if (size(goodConfigs,1)>10000)
                inds = find(distances <= thresh);
                goodConfigs = configs(inds(1:100), :); % all with the same error exactly - probably equivalent
        end
end  

tooHighPercentage = (orig_percentage > 0.05);
veryLowPercentage = (orig_percentage < 0.01);
extremelyHighPercentage = (orig_percentage > 0.2);

if (~isempty(bestGridVec))
    [exists,bestGridInd] = IsMemberApprox(goodConfigs,bestGridVec,1000*eps);
    if (~exists)
        disp('problem with configs');
    end
end
end

function [xs, ys] = getPixelSample(mask, numPoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locs = find(mask);
ind = randi(length(locs), [1,numPoints]);
[ys,xs] = ind2sub(size(mask),locs(ind));
ys = ys';
xs = xs';
end

function [xs, ys] = getPixelSamplePatch(mask, numPoints)
% use 3x3 patch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numCenters = round(numPoints/7); % 补偿一点由于重叠导致的点不够
locs = find(mask);
ind = randi(length(locs), [1,numCenters]);
[ys,xs] = ind2sub(size(mask),locs(ind));
ys = [ys-1; ys-1; ys-1; ys; ys; ys; ys+1; ys+1; ys+1];
xs = [xs-1; xs; xs+1; xs-1; xs; xs+1; xs-1; xs; xs+1];
inliers = (xs>0 & xs<size(mask,2)) & (ys>0 & ys<size(mask,1));
ys = ys(inliers);
xs = xs(inliers);
newInd = sub2ind(size(mask), ys, xs);
newInd = unique(newInd);
inliersInd = mask(newInd)==1;
newInd = newInd(inliersInd);
[ys,xs] = ind2sub(size(mask),locs(newInd));

ys = ys';
xs = xs';
end