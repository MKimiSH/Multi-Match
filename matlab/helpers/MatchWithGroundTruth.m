function [matchInds, corrOverlapErrors] = MatchWithGroundTruth(affines, gts)

nAffs = length(affines);
nGts = length(gts);

matchInds = zeros(length(affines), 1);
isMatched = zeros(length(gts), 1);
overlapErrors = zeros(length(affines), length(gts));


for i=1:nAffs
    for j=1:nGts
        overlapErrors(i,j) = CalculateOverlapError(affines(i), gts(j));
    end
end

corrOverlapErrors = zeros(length(affines),1);
for i = 1:nAffs
    [minOE, idx] = min(overlapErrors(i,:));
    if isMatched(idx)
        corrOverlapErrors(i) = 100000;
        continue;
    else
        corrOverlapErrors(i) = minOE;
        if(minOE < 0.5)
            isMatched(idx) = 1;
            matchInds(i) = idx;
        end
        
    end
end