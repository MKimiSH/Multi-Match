function [configs] = MaskOutlierConfigs(configs, mask)
% Check whether (the translations of) _configs_ are within _mask_

[h, w] = size(mask);

r2x = 0.5*(w-1);
r2y = 0.5*(h-1);

rtx = round(configs(:, 1) + r2x); % tx, ty
rty = round(configs(:, 2) + r2y);
withinArea = (rtx > 0) & (rtx <= w) ...
             & (rty > 0) & (rty <= h); 
         % 不确定这个是不是好的评判标准，先试试
rtx = rtx(withinArea==1);
rty = rty(withinArea==1);
configs = configs(withinArea, :);

indQuery = sub2ind(size(mask), rty, rtx);

withinMask = mask(indQuery) == 1;
% withinMaskInd = find(withinMask);
configs = configs(withinMask, :);

end