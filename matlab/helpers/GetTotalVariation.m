function totalVariation = GetTotalVariation(I1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = size(I1,3);
totalVariation = 0;

for dim = 1 : d
    I1_dim = I1(:,:,dim);
    
    maxVals = ordfilt2(I1_dim,9,true(3));
    minVals = ordfilt2(I1_dim,1,true(3));
    perimeter = max(maxVals-I1_dim,I1_dim-minVals);
    totalVariation = totalVariation + mean2(perimeter);    
end

totalVariation = totalVariation/d;

