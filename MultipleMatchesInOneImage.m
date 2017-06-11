function [] = MultipleMatchesInOneImage(img, tpl, mconfigs, gtconfigs)
% show multiple matches in one image

figure, imshow(img), hold on

[h1,w1,d] = size(tpl);
[h2,w2,d] = size(img);

r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);

for i=1:size(mconfigs, 1)
    cornersX = [1 w1 w1 1];
    cornersY = [1 1 h1 h1];
    transMat = CreateAffineTransformation(mconfigs(i,:));
    a = transMat;
    tx = a(1,3);
    ty = a(2,3);
    
    a2x2 = a(1:2,1:2);
    
    cornersA = a2x2*[cornersX-(r1x+1);cornersY-(r1y+1)];
    cornerAxs = round(cornersA(1,:) + (r2x+1)  + a(1,3));
    cornerAys = round(cornersA(2,:) + (r2y+1)  + a(2,3));
    plot([cornerAxs cornerAxs(1)],[cornerAys cornerAys(1)], '.-g', 'LineWidth', 1.2);
end

for j=1:size(gtconfigs, 1)
    cornersX = [1 w1 w1 1];
    cornersY = [1 1 h1 h1];
    transMat = CreateAffineTransformation(gtconfigs(j,:));
    a = transMat;
    tx = a(1,3);
    ty = a(2,3);
    
    a2x2 = a(1:2,1:2);
    
    cornersA = a2x2*[cornersX-(r1x+1);cornersY-(r1y+1)];
    cornerAxs = round(cornersA(1,:) + (r2x+1)  + a(1,3));
    cornerAys = round(cornersA(2,:) + (r2y+1)  + a(2,3));
    plot([cornerAxs cornerAxs(1)],[cornerAys cornerAys(1)], '.-y', 'LineWidth', 1.2);
end


end