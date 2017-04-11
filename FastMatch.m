function [bestConfig,bestTransMat,sampledError] = FastMatch(template,img,... % mandatory
    templateMask, epsilon,delta,photometricInvariance, searchRange) %optional                                                                        
%FastMatch Find an approximate template match under affine transformations.
%   FastMatch(template,img,... % mandatory
%    [epsilon=0.15],[delta=0.25],[photometricInvariance=false],[searchRange]) %optional   
%  The img and template should both be of class ''double'' (in the range [0,1])')
% searchRange is a struct with fields:
%     'minScale': should be in the range [0,1] default = 0.5 (shrink by factor of at most 2 in x and y scales)
%     'maxScale': should be in the range [1,5) default = 2 (expand by factor of at most 2 in x and y scales)
%     'minRotation': should be in the range [-pi,0] default = -pi (direction is clockwise)
%     'maxRotation': should be in the range [0,pi] default = pi (direction is clockwise)
%     'minTx','maxTx','minTy','maxTy': translations relative to the centers of both images


%% set default values for optional variables
if (~exist('epsilon','var') || isempty(epsilon))
    epsilon = 0.15;
end
if (~exist('delta','var') || isempty(delta))
    delta = 0.25;
end
if (~exist('photometricInvariance','var') || isempty(photometricInvariance))
    photometricInvariance = 0;
end
if (~exist('templateMask','var') || isempty(templateMask))
    templateMask = ones(size(template,1), size(template,2));
end


%% formatting the image and template
img = MakeOdd(img);
template = MakeOdd(template);
templateMask = MakeOdd(templateMask);

if ( ~strcmp(class(img),'double') || ~strcmp(class(template),'double') || ...
        (min([img(:);template(:)])<-0.1) || (max([img(:);template(:)])>1.1) ) %#ok<STISA>
    error('FastMatch: img and template should both be of class ''double'' (in the approx range of [0,1])');
end

if ( size(img,3) ~= size(template,3))
    error('img and template should both be of same dimension per pixel');
end


%% image dimensions
[h1,w1,d1] = size(template);
[h2,w2,d2] = size(img);
r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);


%% search range

if ~exist('searchRange','var')
    searchRange = [];
end
if ~isfield(searchRange,'minScale'), searchRange.minScale = 0.5; end
if ~isfield(searchRange,'maxScale'), searchRange.maxScale = 2; end
if ~isfield(searchRange,'minRotation'), searchRange.minRotation = -pi; end
if ~isfield(searchRange,'maxRotation'), searchRange.maxRotation = pi; end
if ~isfield(searchRange,'minTx'), searchRange.minTx = -(r2x-r1x*searchRange.minScale); end
if ~isfield(searchRange,'maxTx'), searchRange.maxTx = r2x-r1x*searchRange.minScale; end
if ~isfield(searchRange,'minTy'), searchRange.minTy = -(r2y-r1y*searchRange.minScale); end
if ~isfield(searchRange,'maxTy'), searchRange.maxTy = r2y-r1y*searchRange.minScale; end

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
minTx = max(searchRange.minTx,-(r2x-r1x*minScale));
maxTx = min(searchRange.maxTx,r2x-r1x*minScale);
minTy = max(searchRange.minTy,-(r2y-r1y*minScale));
maxTy = min(searchRange.maxTy,r2y-r1y*minScale);


%% parametrize the initial grid
[bounds,steps] = GenerateGrid(w1,h1,delta,minTx,maxTx,minTy,maxTy,minRotation,maxRotation,minScale,maxScale);


%% run the actual search!
[bestConfig,bestTransMat,sampledError] = ...
    FindBestTransformation(template,img,bounds,steps,epsilon,delta,photometricInvariance, templateMask);

%%
return

