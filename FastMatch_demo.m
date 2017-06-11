function FastMatch_demo
%%%%%%%%%%%%%%%%%%%%%%%

clc
% clear all
close all
dbstop if error

% adding 2 subdirectories to Matlab PATH
AddPaths

% compiling the relevant Mex-files
CompileMex


%% Select examples to run

example_1 = 0; % "predefined template in a predefined image"
example_2 = 0; % "predefined template with given MASK in a predefined image"
example_3 = 0; % "RANDOM template in a user-selected image"
example_4 = 1; % "USER-SELECTED template in a user-selected image"
example_5 = 0; % "choose a pair of RELATED images and select a template in the first"
example_6 = 0; % "random template in a user-selected image with a limited SEARCH RANGE"

%%

%% example 1: "predefined template in a predefined image"
if example_1
    disp('example 1: Locating a predefined template in an image (thai food)');
    fprintf('======================================================\n');
        
    % reading image and template
    img =       imread([pwd '/example/image.png']);
    img =       im2double(img);
    img = MakeOdd(img);
    template =  imread([pwd './example/template.png']);
    template =  im2double(template);
    template = MakeOdd(template);
        
    % reading the 'ground-truth' transformation (not mandatory)
    load([pwd '/example/OptMat.mat'], 'optMat');
    
    % FastMatch run
    [bestConfig,bestTransMat,sampledError] = FastMatch(template,img);
        
    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,optMat,'example 1');
    
    fprintf('example 1 - optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);
    fprintf('example 1: finished\n\n');
end


%% example 2: "predefined template with given MASK in a predefined image"
if example_2
    disp('example 2: Locating a predefined MASKED template in an image (thai food)');
    fprintf('======================================================\n');
        
    % reading image and template
    img =       imread([pwd '/example/image.png']);
    img =       im2double(img);
    img = MakeOdd(img);

    template =  imread([pwd './example/template.png']);
    template =  im2double(template);
    template = MakeOdd(template);
    
    % Demo for template mask
    s1 = size(template,1);
    s2 = size(template,2);
    templateMaskY = 1:s1;
    templateMaskX = 1:s2;
    templateMask = templateMaskY' * templateMaskX;
    templateMask = (templateMask > sqrt(s1)*sqrt(s2)) & (templateMask < 6*sqrt(s1)*sqrt(s2));
    
    % reading the 'ground-truth' transformation (not mandatory)
    load([pwd '/example/OptMat.mat'], 'optMat');
    
    % FastMatch run
    [bestConfig,bestTransMat,sampledError] = FastMatch(template,img,templateMask);
    
    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,optMat,'example 2', templateMask);
    
    fprintf('example 2 -  optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);    
    fprintf('example 2: finished\n\n');
       
end



%% example 3: "RANDOM template in a user-selected image"
if example_3    
    
    disp('example 3: locating a randomly generated template in a user-selected image');
    fprintf('======================================================\n');    
    
    % reading image and template
    img = SelectAnImage('Choose an image for example 3','imgs\*.*');
    img = MakeOdd(img);

    [h,w,d] = size(img);
    minDim2 = min(h,w); % minimal dimension
    sizeFact = 2.5; % maximal relative size of template dimension w.r.t. image dimension
    n1 = ceil(minDim2/sizeFact);
    
    [template,optMat] = GenerateRandomAffineTemplate(img,n1);

    % FastMatch run
    [bestConfig,bestTransMat,sampledError] = FastMatch(template,img);
    
    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,optMat,'example 3');
    
    fprintf('example 3 -  optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);
    fprintf('example 3: finished\n\n');  
       
end



%% example 4: "USER-SELECTED template in a selected image"
if example_4
        
    disp('example 4: locating a user selected template in an image');
    fprintf('======================================================\n');
    
    % reading image and template
    img = SelectAnImage('Choose an image for example 4','imgs\*.*');
    img = MakeOdd(img);

    fprintf('You will need to mark 3 corners of a parallelogram using a triangle\n');
    fprintf('(don''t make smaller than original)\n\n');
    fprintf('(move/stretch/rotate the suggested triangle and DOUBLE-CLICK it)\n\n');
    [template,optMat] = GenerateUserSelectedTemplate(img,'Example ');
        
    % FastMatch run
    [bestConfig,bestTransMat,sampledError] = FastMatch(template,img);
    
    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,optMat,'example 4');
    
    fprintf('example 4 -  optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);   
    fprintf('example 4: finished\n\n');
    
end


%% example 5: "Choose a pair of RELATED images and select a template in the first"
if example_5
    disp('example 5: "Choose related images" - locating a user selected template from one image in another related image');
    fprintf('======================================================\n');
    
    disp('Loading a target image and a template image for example 5...');
    fprintf('(make sure to choose target and template with same index (1-5))\n\n');
    fprintf('(then - move or enlarge the suggested template and DOUBLE-CLICK it)\n\n');
    
    % reading image and template    
    img = SelectAnImage('Choose an image to match with the template','imgPairs\*.*');
    img = MakeOdd(img);

    template_img = SelectAnImage('Choose an image to match with the template','imgPairs\*.*');
    template_img = MakeOdd(template_img);

    [template] = GenerateUserSelectedTemplateForImagePair(template_img,img);
    
%     ShowInstance(template,img,'example 5');
    
    % FastMatch run
   [bestConfig,bestTransMat,sampledError] = FastMatch(template,img);
    
    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,[],'example 5');
    % Note that there's no "Ground truth"
    
    fprintf('example 5 -  optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);   
    fprintf('example 5: finished\n\n');
    
end


%% example 6: "Random template in a user-selected image with a limited SEARCH RANGE"
if example_6    
    
    disp('example 6: locating a randomly generated template in a user-selected image with a limited SEARCH RANGE');
    fprintf('======================================================\n');
    
    
    % search range (Note: you can specify a subset of the fields - the rest are assigned defaults)
    searchRange.minScale = 0.4;
    searchRange.maxScale = 2.4;
    searchRange.minRotation = -pi/16;
    searchRange.maxRotation = pi/2;
    searchRange.minTx = -50;
    searchRange.maxTy = 150;
    
    
    % reading image and template
    img = SelectAnImage('Choose an image for example 6','imgs\*.*');
    img = MakeOdd(img);

    [h,w,d] = size(img);
    minDim2 = min(h,w); % minimal dimension
    sizeFact = 3; % relative size of template dimension w.r.t. image dimension
    n1 = ceil(minDim2/sizeFact);
    % Note: random template is generated within the search range
    [template,optMat] = GenerateRandomAffineTemplate(img,n1,searchRange);
    

    % FastMatch run
    [bestConfig,bestTransMat,sampledError] = FastMatch(template,img,[],[],[],[],searchRange);

    % Visualize result
    [optError,fullError,overlapError] = MatchingResult(template,img,bestTransMat,optMat,'example 6');
 
    fprintf('example 6 -  optError: %.4f (%.2f GLs), fullError: %.4f (%.2f GLs), overlapError: %.1f%%\n',...
        optError,256*optError,fullError,256*fullError,100*overlapError);  
    fprintf('example 6: finished\n\n');
    
end


%%

fprintf('END OF DEMO!\n\n');

keyboard

return
end



% % % % % % % Helper functions % % % % % % % %

function [template,optMat] = GenerateRandomAffineTemplate(img,n1,searchRange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pick a random transformation and continue if it is in bounds

[h,w,d] = size(img);

if (mod(n1,2)==0)
    n1 = n1 - 1;
end

r1 = 0.5*(n1-1);
rx = 0.5*(w-1);
ry = 0.5*(h-1);

%% affine search limits

if ~exist('searchRange','var')
    searchRange = [];
end
if ~isfield(searchRange,'minScale'), searchRange.minScale = 0.5; end
if ~isfield(searchRange,'maxScale'), searchRange.maxScale = 2; end
if ~isfield(searchRange,'minRotation'), searchRange.minRotation = -pi; end
if ~isfield(searchRange,'maxRotation'), searchRange.maxRotation = pi; end
if ~isfield(searchRange,'minTx'), searchRange.minTx = -(rx-r1*searchRange.minScale); end
if ~isfield(searchRange,'maxTx'), searchRange.maxTx = rx-r1*searchRange.minScale; end
if ~isfield(searchRange,'minTy'), searchRange.minTy = -(ry-r1*searchRange.minScale); end
if ~isfield(searchRange,'maxTy'), searchRange.maxTy = ry-r1*searchRange.minScale; end

% check ranges
assert(searchRange.minScale >=0 && searchRange.minScale <=1);
assert(searchRange.maxScale >=1 && searchRange.maxScale <=5);
assert(searchRange.minRotation >=-pi && searchRange.minRotation <=0);
assert(searchRange.maxRotation >=0 && searchRange.maxRotation <=pi);

% copy params
minScale = searchRange.minScale;
maxScale = searchRange.maxScale;
minRotation = searchRange.minRotation;
maxRotation = searchRange.maxRotation;
minTx = max(searchRange.minTx,-(rx-r1*minScale));
maxTx = min(searchRange.maxTx,rx-r1*minScale);
minTy = max(searchRange.minTy,-(ry-r1*minScale));
maxTy = min(searchRange.maxTy,ry-r1*minScale);

rangeTx = maxTx - minTx;
rangeTy = maxTy - minTy;
rangeR = 2*pi;
rangeS = maxScale - minScale;

attempt = 0;
while(1)
    attempt = attempt + 1;
    randVec = rand(1,6);
    randAffVec = [minTx,minTy,-pi,minScale,minScale,-pi] + randVec.*[rangeTx,rangeTy,rangeR/4,rangeS,rangeS,rangeR];
    
    % check rotation range:
    % total rotation in the range [0,2*pi]
    totalRot = mod(randAffVec(3)+randAffVec(6),2*pi);
    % total rotation in the range [-pi,pi]
    totalRot(totalRot>pi) = totalRot(totalRot>pi) - 2*pi;
    % filtering
    inRange = totalRot>=minRotation && totalRot<=maxRotation;
    if ~inRange
        continue
    end

    %             randAffVec = [0 0 0 1 1 0];
    randAffMat = CreateAffineTransformation(randAffVec);
    %         randAff = maketform('affine',randAffMat');
    cornersX = [1 n1 n1 1];
    cornersY = [1 1 n1 n1];
    randAff2x2 = randAffMat(1:2,1:2);
    cornersTest = (randAff2x2*[cornersX-(r1+1);cornersY-(r1+1)]);
    cornersTestxs = round(cornersTest(1,:) + (rx+1)  + randAffMat(1,3));
    cornersTestys = round(cornersTest(2,:) + (ry+1)  + randAffMat(2,3));
    
    % check if mapped into I2 bounds
    if (isequal(BoundBy(cornersTestxs,1,w),cornersTestxs) && (isequal(BoundBy(cornersTestys,1,h),cornersTestys)))
        %% Construct the template I1
        % get the inverse transform
        centeredTL3 = [1,1; n1,1; n1,n1] - (r1+1); % BL,TL,TR
        p = [cornersTestxs' cornersTestys'];
        aff = cp2tform(centeredTL3,p(1:3,:),'affine');
        
        template = imtransform(img,fliptform(aff),'nearest','XData',[1 n1] - (r1+1),'YData',[1 n1] - (r1+1));
        % verify that it isn't just flat
        if (std(double(template(:)))>0.1)
            break;
        end
    end
end

optMat = randAffMat;

return;
end

function [template,optMat] = GenerateUserSelectedTemplate(img,prefixName,notQuiet)

[h,w] = size(img);

minDim2 = min(h,w); % minimal dimension

sizeFact = 3; % relative size of template dimension w.r.t. image dimension

n1 = ceil(minDim2/sizeFact);

n12 = round(minDim2/2);
n34 = n12+n1-1;


r1 = 0.5*(n1-1);
rx = 0.5*(w-1);
ry = 0.5*(h-1);


%% get parallelogram in I2 from the user and transform it to a square (= query image I1)
tempfig = figure; imshow(img);
title('adjust (dont enlarge) the triangle (half parallelogram) and them DOUBLE CLICK it');
set(gcf,'name',[prefixName ': template and image']);

hand = impoly(gca,[n34,n12; n12,n12; n12,n34]);
p = wait(hand);

p1 = [p(1,1),p(1,2)];
p2 = [p(2,1),p(2,2)];
p3 = [p(3,1),p(3,2)];
mid = 0.5*(p1+p3);
opp = p2 + 2*(mid-p2);



% Sides of the parallelogream
s1 = norm(p1-p2);
s2 = norm(p3-p2);
smallestRelation = min(s1/n1, s2/n1);
largestRelation  = max(s1/n1, s2/n1);
if ((smallestRelation < 0.5) || (largestRelation > 2))
    warning('Template chosen has dimesions too small or too large for demo parameters. You can change the value of the variable "sizeFact" in the FastMAtch_demo file.'); %#ok<WNTAG>
end
hold off;


centeredTL3 = [n1,1; 1,1; 1,n1] - (r1+1); % BL,TL,TR
aff = cp2tform(centeredTL3,p,'affine');

template = imtransform(img,fliptform(aff),'bicubic','XData',[1 n1] - (r1+1),'YData',[1 n1] - (r1+1));

%% extract the parameters of the transformation (for ground-truth visualization)
centeredP = [p(:,1)-(rx+1), p(:,2)-(ry+1)];
aff = cp2tform(centeredTL3,centeredP,'affine');

% aff.tdata.T = [1 0 0; 0 1 0; 0 0 1];
optMat = aff.tdata.T';
centerPoint = [r1+1;r1+1];
tcenterpoint = optMat(1:2,1:2)*centerPoint;
vecx = [-rx,rx];
vecy = [-ry,ry];

I1MappedByOPT = imtransform(template,aff,'bicubic','xdata',tcenterpoint(1)+vecx,'ydata',tcenterpoint(2)+vecy,'size',size(img));

% show the images
close(tempfig);
if exist('notQuiet','var') && notQuiet
    figure;
    set(gcf,'name',[prefixName ': template and image']);
    subplot 121; imshow(template); title('template');
    subplot 122; imshow(img); title('desired location in target img');
    
    subplot 122; hold on; plot([p1(1),p2(1),p3(1),opp(1),p1(1)],[p1(2),p2(2),p3(2),opp(2),p1(2)],'*-g');
end
return;
end



function [template] = GenerateUserSelectedTemplateForImagePair(I1,I2)

%% display images
[h1,w1,d1] = size(I1);
hf = figure
subplot 121; axis on;   imshow(I1);
subplot 122; axis on;   imshow(I2);

subplot 121; hold on;
h = imrect(gca,[w1/2 h1/2 w1/6 w1/6]);
p = wait(h); % p(4) = p(3);
p = round(p);
tempXmin = p(1);
tempYmin = p(2);
tempW = p(3);
tempH = p(4);

%% define template
template = I1(tempYmin:tempYmin+tempH-1,tempXmin:tempXmin+tempW-1);
close(hf);
return;
end


function img = SelectAnImage(str,pref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imgFileName,PathName,~] = uigetfile('*','Choose an image to match with the template',pref);
FullName = [PathName,imgFileName];

img = imread(FullName);

[~,~,d2] = size(img);
if (d2>1)
    img = rgb2gray(img);
end
img = im2double(img);
return
end
