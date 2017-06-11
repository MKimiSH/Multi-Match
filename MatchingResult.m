function [optError,fullError,overlapError, figHandle] = MatchingResult(I1,I2,a,optA,prefixName,templateMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('templateMask','var'))
    templateMask = ones(size(I1,1), size(I1,2));
end

I1(templateMask==0) = 0;
isGrayscale = (size(I1,3)==1);

tx = a(1,3);
ty = a(2,3);

[h1,w1,d] = size(I1);
[h2,w2,d] = size(I2);

r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);


% take the complete set of pixels (not a sample)
% [xs,ys] = meshgrid(1:w1,1:h1);
[ys,xs]=find(templateMask);
xs = xs(:)'; 
ys = ys(:)';


% mapping by A
a2x2 = a(1:2,1:2);
targetPoints = a2x2*[xs-(r1x+1);ys-(r1y+1)];
txs = round(targetPoints(1,:) + (r2x+1)  + tx);
tys = round(targetPoints(2,:) + (r2y+1)  + ty);
insideInds = find(txs>0 & txs<(w2+1) & tys>0 & tys<(h2+1) );


% mapping by OPT
if (~isempty(optA))
    a2x2_OPT = optA(1:2,1:2);
    targetPoints_OPT = (a2x2_OPT*[xs-(r1x+1);ys-(r1y+1)]);
    tx_OPT = optA(1,3);
    ty_OPT = optA(2,3);
    txs_OPT = round(targetPoints_OPT(1,:) + (r2x+1)  + tx_OPT);
    tys_OPT = round(targetPoints_OPT(2,:) + (r2y+1)  + ty_OPT);
    insideInds_OPT = find(txs_OPT>0 & txs_OPT<(w2+1) & tys_OPT>0 & tys_OPT<(h2+1) );
end

cornersX = [1 w1 w1 1];
cornersY = [1 1 h1 h1];
if (~isempty(optA))
    optA2x2 = optA(1:2,1:2);
    cornersOPT = optA2x2*[cornersX-(r1x+1);cornersY-(r1y+1)];
    cornerOPTxs = round(cornersOPT(1,:) + (r2x+1)  + optA(1,3));
    cornerOPTys = round(cornersOPT(2,:) + (r2y+1)  + optA(2,3));
end
cornersA = a2x2*[cornersX-(r1x+1);cornersY-(r1y+1)];
cornerAxs = round(cornersA(1,:) + (r2x+1)  + a(1,3));
cornerAys = round(cornersA(2,:) + (r2y+1)  + a(2,3));

xLimits = BoundBy([min(cornerAxs),max(cornerAxs)],1,w2);
yLimits = BoundBy([min(cornerAys),max(cornerAys)],1,h2);

if (isempty(insideInds))
    error('in CompareToOPTandSaveWorkspace - THIS SHOULDN''T HAPPEN!!!');
end

%% full (l_1 distance) error
numPoints = length(xs);
sourceInds = sub2ind([h1,w1],ys(insideInds),xs(insideInds));
targetInds = sub2ind([h2,w2],tys(insideInds),txs(insideInds));

badMatches =  sum(abs(I1(sourceInds) - I2(targetInds)));
fullError = (badMatches + numPoints - length(insideInds))/numPoints;



if (~isempty(optA))
    %% full (l_1 distance) error for OPT
    numPoints = length(xs);
    sourceInds_OPT = sub2ind([h1,w1],ys(insideInds_OPT),xs(insideInds_OPT));
    targetInds_OPT = sub2ind([h2,w2],tys_OPT(insideInds_OPT),txs_OPT(insideInds_OPT));

    badMatches_OPT =  sum(abs(I1(sourceInds_OPT) - I2(targetInds_OPT)));
    optError = (badMatches_OPT + numPoints - length(insideInds_OPT))/numPoints;

    %% overlap error
    try
    [cornerOPTxs,cornerOPTys] = poly2cw(cornerOPTxs,cornerOPTys);
    [cornerAxs,cornerAys] = poly2cw(cornerAxs,cornerAys);
    [unionXs, unionYs] = polybool('union', cornerOPTxs,cornerOPTys,cornerAxs,cornerAys);
    [interXs, interYs] = polybool('intersection',  cornerOPTxs,cornerOPTys,cornerAxs,cornerAys);

    areaUnion = polyarea(unionXs, unionYs);
    areaIntersect = polyarea(interXs, interYs);
    catch
        warning('Probably missing ''map'' toolbox - Not calculating intersection score!!!')
        areaUnion = -999;
        areaIntersect = -999;
    end

    if (areaIntersect==0)
        overlapError = 1;
    else
        overlapError = 1 - (areaIntersect/areaUnion);
    end
else
    optError = -1;
    overlapError = -1;
end

%% an alternative way of computing the full error
if isGrayscale
    FM_mexDist = EvaluateConfigs_mex(I1',I2',a([1 4 7 2 5 8])',int32(xs),int32(ys),0);
    if (~isempty(optA))
        OPT_mexDist = EvaluateConfigs_mex(I1',I2',optA([1 4 7 2 5 8])',int32(xs),int32(ys),0);
    else
        OPT_mexDist = -1;
    end
else
    FM_mexDist = EvaluateConfigsVectorized_mex(permute(I1,[3,2,1]),permute(I2,[3,2,1]),a([1 4 7 2 5 8])',int32(xs),int32(ys),0);
    if (~isempty(optA))
        OPT_mexDist = EvaluateConfigsVectorized_mex(permute(I1,[3,2,1]),permute(I2,[3,2,1]),optA([1 4 7 2 5 8])',int32(xs),int32(ys),0);
    else
        OPT_mexDist = -1;
    end
end
fprintf('(Alternative Mex error computations) OPT_mexDist: %.4f, FM_mexDist: %.4f\n',OPT_mexDist,FM_mexDist);
    
%% transform I1 using A
tran = maketform('affine',a');
centerpoint = [r1x+1;r1y+1];
tcenterpoint = a(1:2,1:2)*centerpoint;
I2rangeX = [-r2x,r2x];
I2rangeY = [-r2y,r2y];

% transformed I1
transtarget = imtransform(I1,tran,'bilinear','xdata',tcenterpoint(1)+I2rangeX,'ydata',round(tcenterpoint(2)+I2rangeY),'size',[size(I2,1),size(I2,2)]);

% a mask on ones
%I3 = I1;
%I3(:) = 1;
I3 = templateMask;
% transformed mask
transI3 = imtransform(I3,tran,'bilinear','xdata',tcenterpoint(1)+I2rangeX,'ydata',round(tcenterpoint(2)+I2rangeY),'size',[size(I2,1),size(I2,2)]);
transI3 = repmat(transI3,[1,1,d]);
transI3mask = transI3 > 0.5;


%% summarizing figure (image)
tic
fullscreen = get(0,'ScreenSize');
figHandle = figure();
set(gcf,'Position',[0.15*fullscreen(3) 0.1*fullscreen(4) 0.8*fullscreen(3) 0.75*fullscreen(4)]);
% GL.hh = hh;
set(gcf,'color','w');
set(gcf,'name',[prefixName, ': FM result']);

%I1
subplot 221; hold off;
% I1padded = padarray(I1,[max(0,h2-h1),max(0,floor(0.5*(w2-w1)))],1,'post'); 
% I1padded = padarray(I1padded,[0,max(0,ceil(0.5*(w2-w1)))],1,'pre'); 
white_I1 = I1;
white_I1(~templateMask) = 1;
I1padded_the_size_of_I2 = padarray(white_I1,[round((h2-h1)/2),round((w2-w1)/2)],1);
imshow(I1padded_the_size_of_I2); hold on;
title(['TemplateSize: ' num2str(h1) 'x' num2str(w1) ...
        '.   TV: ' num2str(GetTotalVariation(I1),'%.3f') '.  STD: ' num2str(std(I1(:)),'%.3f')]);

% I2
subplot 222; hold off;
imshow(I2); hold on;
secondSTR = sprintf('\nwarped template (magenta)');
if (~isempty(optA))
    plot([cornerOPTxs cornerOPTxs(1)],[cornerOPTys cornerOPTys(1)], '*-g');
    secondSTR = sprintf('%s, ground truth (green)',secondSTR);
end
plot([cornerAxs cornerAxs(1)],[cornerAys cornerAys(1)], '*-m');
% title(['TargetSize: ' num2str(h2) 'x' num2str(w2) '.     outside: ' num2str((numPoints-length(insideInds))/numPoints,'%.3f') ' L1-dist: ' num2str(fullError,'%.4f') ' overlapErr: ' num2str(overlapError,'%.4f')]);
title(['TargetSize: ' num2str(h2) 'x' num2str(w2) '.   L1-dist: ' num2str(fullError,'%.3f') ...
    '   overlapErr: ' num2str(overlapError,'%.3f') secondSTR]);


% query result
subplot 223; hold off;
maskedtranstarget = ~transI3mask.*ones(size(transI3mask)) + transI3mask.*transtarget;
imshow(maskedtranstarget(yLimits(1):yLimits(end),xLimits(1):xLimits(end),:)); hold on;
title('the warped template');

% ground truth
subplot 224; hold off;
transI3mask = double(transI3mask);
target =  ~transI3mask.*ones(size(transI3mask)) + transI3mask.*I2;

imshow(target(yLimits(1):yLimits(end),xLimits(1):xLimits(end), :)); hold on;
title('background target');


return







