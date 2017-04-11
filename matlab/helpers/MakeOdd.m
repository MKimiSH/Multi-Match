function [I,croppedH,croppedW] = MakeOdd(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k: make the image dims odd!

[h,w,d] = size(I);

    croppedH = 0;
    croppedW = 0;

if (mod(h,2)==0)
    croppedH = 1;
    I = I(1:end-1,:,:);
end

if (mod(w,2)==0)
    croppedW = 1;
    I = I(:,1:end-1,:);
end

if (mod(d,2)==0)
    I = I(:,:,1:end-1);
end