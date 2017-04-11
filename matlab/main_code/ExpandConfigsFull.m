function expandedConfigs = ExpandConfigsFull(configs,steps,level,deltaFact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fact = deltaFact^level;
halfstep_tx = steps.tx/fact;
halfstep_ty = steps.ty/fact;
halfstep_r = steps.r/fact;
halfstep_s = steps.s/fact;

npoints = 3^6;
additions = zeros(npoints,6);
i = 0;
for tx = -halfstep_tx : halfstep_tx : halfstep_tx
    for ty = -halfstep_ty : halfstep_ty : halfstep_ty
        for r2 = -halfstep_r : halfstep_r : halfstep_r
            for sx = -halfstep_s : halfstep_s : halfstep_s
                for sy = -halfstep_s : halfstep_s : halfstep_s
                    for r1 = -halfstep_r : halfstep_r : halfstep_r
                        i = i + 1;
                        additions(i,:) = [tx,ty,r2,sx,sy,r1];
                    end
                end
            end
        end
    end
end

numConfigs = size(configs,1);
expanded= repmat(configs,[npoints,1]);

expandedConfigs = expanded + repmat(additions,[numConfigs,1]);

disp('done');



