function expandedConfigs = ExpandConfigsRandom(configs,steps,level,npoints,deltaFact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fact = deltaFact^level;
halfstep_tx = steps.tx/fact;
halfstep_ty = steps.ty/fact;
halfstep_r = steps.r/fact;
halfstep_s = steps.s/fact;


numConfigs = size(configs,1);
% randvec = 2*(round(rand(npoints*numConfigs,6))-0.5); % random vectors in {-1,1}^6
randvec = floor(3*rand(npoints*numConfigs,6)-1); % random vectors in {-1,0,1}^6
expanded= repmat(configs,[npoints,1]);
ranges = [halfstep_tx,halfstep_ty,halfstep_r,halfstep_s,halfstep_s,halfstep_r];

expandedConfigs = expanded + randvec.*repmat(ranges,[npoints*numConfigs,1]);

% disp('done');



