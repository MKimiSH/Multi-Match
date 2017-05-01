function [configs] = BoundConfigsTrSc(configs, bounds)
% Bound configs inside _bounds_
% to avoid too small scaling...
% No need to bound rotation because it is bounded in FindOneMatch()

btx = configs(:,1) >= bounds.tx(1) & configs(:,1) <= bounds.tx(2);
bty = configs(:,2) >= bounds.ty(1) & configs(:,2) <= bounds.ty(2);
bsx = configs(:,4) >= bounds.s(1) & configs(:,4) <= bounds.s(2);
bsy = configs(:,5) >= bounds.s(1) & configs(:,5) <= bounds.s(2);

boundInds = btx & bty & bsx & bsy;
configs = configs(boundInds, :);

end