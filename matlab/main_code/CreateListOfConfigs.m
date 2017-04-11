function [configs,gridSize] = CreateListOfConfigs(bounds,steps,onlyGridSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx_steps = bounds.tx(1) : steps.tx : (bounds.tx(2) + 0.5*steps.tx); % 'pad' at end of range with an extra sample
ty_steps = bounds.ty(1) : steps.ty : (bounds.ty(2) + 0.5*steps.ty); % 'pad' at end of range with an extra sample
% Rotations ignore the user selected range here - it is handles in FindBestTransformation
r_steps = -pi : steps.r : pi; % no padding since it is a cyclic range
s_steps = bounds.s(1) : steps.s : bounds.s(2) + 0.5*steps.s; % 'pad' at end of range with an extra sample
if (steps.s == 0)
    s_steps = bounds.s(1);
end

% number of steps
ntx_steps = length(tx_steps);
nty_steps = length(ty_steps);
ns_steps = length(s_steps);
nr_steps = length(r_steps);

% second rotation is a special case (can be limited to a single quartile)
quartile1_r_steps = r_steps(r_steps < -pi/2 + steps.r/2);
NR2_steps = length(quartile1_r_steps);

% gridsize
gridSize = ntx_steps*nty_steps*(ns_steps^2)*(nr_steps*NR2_steps);

if (exist('onlyGridSize','var') && onlyGridSize)
    configs = [];
    return
end

% tic
% configs = CreateList(nt_steps,nr_steps,NR2_steps,ns_steps,tx_steps,r_steps,s_steps,gridSize);
% toc
configs_mex = CreateList_mex(ntx_steps,nty_steps,nr_steps,NR2_steps,ns_steps,tx_steps,ty_steps,r_steps,s_steps,int32(gridSize));
configs = configs_mex';


% backup Matlab version
function configs = CreateList(ntx_steps,nty_steps,nr_steps,NR2_steps,ns_steps,tx_steps,ty_steps,r_steps,s_steps,gridSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
configs = zeros(gridSize,6);
gridInd = 0;
for tx_ind = 1 : ntx_steps 
    tx = tx_steps(tx_ind);
%     fprintf('step %d out of %d\n',tx_ind,nt_steps);
    for ty_ind = 1 : nty_steps 
        ty = ty_steps(ty_ind);
        for r1_ind = 1 : nr_steps
            r1 = r_steps(r1_ind);
            for r2_ind =  1 : NR2_steps
                r2 = r_steps(r2_ind);
                for sx_ind =  1 : ns_steps
                    sx = s_steps(sx_ind);
                    for sy_ind =  1 : ns_steps 
                        sy = s_steps(sy_ind);
                        
                        gridInd = gridInd + 1;
                        
                        configs(gridInd,:) = [tx,ty,r2,sx,sy,r1];
                    end
                end
            end
        end
    end
end

