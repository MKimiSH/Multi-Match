function [dbyt, dbys, dbyr] = ShowDistances(distances, configs, bounds, steps)

tx_steps = bounds.tx(1) : steps.tx : (bounds.tx(2) + 0.5*steps.tx); % 'pad' at end of range with an extra sample
ty_steps = bounds.ty(1) : steps.ty : (bounds.ty(2) + 0.5*steps.ty); % 'pad' at end of range with an extra sample
ntx_steps = length(tx_steps);
nty_steps = length(ty_steps);

r_steps = -pi : steps.r : pi; % no padding since it is a cyclic range
s_steps = bounds.s(1) : steps.s : bounds.s(2) + 0.5*steps.s; % 'pad' at end of range with an extra sample

ns_steps = length(s_steps);
nr_steps = length(r_steps);

dbyt = zeros(nty_steps, ntx_steps);
dbys = zeros(ns_steps, ns_steps);
dbyr = zeros(nr_steps, nr_steps);

for sy = 1:ns_steps
    for sx = 1:ns_steps
        dbys(sy,sx) = mean(distances(configs(:,5)==s_steps(sy) & configs(:,4)==s_steps(sx)));
    end
end

for i=1:nty_steps
    for j=1:ntx_steps
        dbyt(i,j) = mean(distances(configs(:,2)==ty_steps(i) & configs(:,1)==tx_steps(j)));
    end
end

for i=1:nr_steps
    for j=1:nr_steps
        dbyr(i,j) = mean(distances(configs(:,3)==r_steps(i) & configs(:,6)==r_steps(j)));
    end
end

figure,
subplot(221), imagesc(dbyt), title('tx-ty');
subplot(222), imagesc(dbys), title('sx-sy');
subplot(223), imagesc(dbyr), title('r1-r2');
end