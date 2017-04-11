function thresh = GetThreshPerDelta(delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on experimentally drawn % p = 0.1341    0.0278

p = [0.1341, 0.0278];
safety = 0.02;
thresh = p(1)*delta+p(2)-safety; %+0.01;