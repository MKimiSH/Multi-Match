function [overlapError] = CalculateOverlapError(a,optA)

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



if (~isempty(optA))
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
    overlapError = 100000; % I don't want to use -1!!!!!!!
end

end