function [img, tpl, factor] = AdjustSize(img, tpl, pz)
% adjust the sizes of image and template
% factor is the difference between the resize factors of _img_ and _tpl_.
% scale range should multiply _factor_
% _img_ -> numel(img) < 1500000(~300x400)
% _tpl_ -> numel(tpl) >= 1000(~30x30)

% imgFactor = 1;
% tplFactor = 1;
factor = 1;

[ih, iw, ~] = size(img);
[th, tw, ~] = size(tpl);

imgFactor = 1450000/(ih*iw);
if(th*tw*imgFactor >= 1000)
    tplFactor = imgFactor;
else
    tplFactor = 1000/(th*tw);
    factor = imgFactor/tplFactor; % scale range should multiply _factor_
end

img = imresize(img, imgFactor);
tpl = imresize(tpl, tplFactor);

img = MakeDivided(img, pz);
tpl = MakeDivided(tpl, pz);

end

function [im] = MakeDivided(im, pz)

sz = size(im);
mt = mod(sz(1:2), pz);

if(any(mt))
    im = im(1:end-mt(1), 1:end-mt(2), :);
end

end