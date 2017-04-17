function [mask, configs] = MaskNewMatch(mask, tpl, configs, aff)
% Delete the points in _mask_ and _configs_ that belong to the transformed _tpl_ by _aff_

tmask = ones(size(tpl,1), size(tpl,2)); % 也许可以不全mask掉，可以只mask一部分（如中心部分）。
[h, w] = size(mask);
[th, tw] = size(tmask);
r2x = 0.5*(w-1);
r2y = 0.5*(h-1);
r1x = 0.5*(tw-1);
r1y = 0.5*(th-1);

[tmaskTrans, xd, yd] = imtransform(tmask, aff);
xrange = round(xd + r2x - r1x);
yrange = round(yd + r2y - r1y);

% border
if(xrange(1)<=0) 
    tmaskTrans = tmaskTrans(:, 2-xrange(1):end);
    xrange(1) = 1;
end
if(xrange(2)>w)
    tmaskTrans = tmaskTrans(:, 1:end-(xrange(2)-w));
    xrange(2) = w;
end
if(yrange(1)<=0)
    tmaskTrans = tmaskTrans(2-yrange(1):end, :);
    yrange(1) = 1;
end
if(yrange(2)>h)
    tmaskTrans = tmaskTrans(1:end-(yrange(2)-h), :);
    yrange(2) = h;
end

tmaskTrans = tmaskTrans<0.5;
szTtr = size(tmaskTrans);

maskPart = mask(yrange(1):yrange(1)+szTtr(1)-1, xrange(1):xrange(1)+szTtr(2)-1);
maskPart = maskPart & tmaskTrans;
mask(yrange(1):yrange(1)+szTtr(1)-1, xrange(1):xrange(1)+szTtr(2)-1) = maskPart;

configs = MaskOutlierConfigs(configs, mask);

end