function counter = GetGrid(sourceW,sourceH,delta,c,targetW,targetH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridStepW = sourceW*delta;
gridStepH = sourceH*delta;
xgrid = (-targetW/2):gridStepW:targetW/2;
ygrid = (-targetH/2):gridStepH:targetH/2;

I1_area = sourceW*sourceH;
extended_I1_area = (sourceW+gridStepW)*(sourceH+gridStepH);
reduced_I1_area = (sourceW-gridStepW)*(sourceH-gridStepH);


p1xIdx = 0;
counter = 0;
opcounter = 0;
for p1x = -targetW/2 : gridStepW : targetW/2
    p1xIdx = p1xIdx + 1;
    fprintf('x index %d out of %d\n',p1xIdx,length(xgrid));
    for p1y = -targetH/2 : gridStepH : targetH/2
        for p2x = max(-targetW/2,p1x-sourceW*c^2) : gridStepW : min(targetW/2,p1x+sourceW*c^2)
            for p2y = max(-targetH/2,p1y-sourceH*c^2) : gridStepH : min(targetH/2,p1y+sourceH*c^2)
                for p3x = max(-targetW/2,p1x-sourceW*c^2) : gridStepW : min(targetW/2,p1x+sourceW*c^2)
                    for p3y = max(-targetH/2,p1y-sourceH*c^2) : gridStepH : min(targetH/2,p1y+sourceH*c^2)
                        p2xShift = p2x-p1x;
                        p3yShift = p3y-p1y;
                        p2yShift = p2y-p1y;
                        p3xShift = p3x-p1x;
                        para_area = abs(p2xShift*p3yShift-p2yShift*p3xShift);
                        if ((para_area > reduced_I1_area/c^2)&& (para_area < extended_I1_area*c^2))
                            counter = counter + 1;
                        else
                            opcounter = opcounter + 1;
                        end
                    end
                end
            end
        end
    end
end

counter
opcounter


