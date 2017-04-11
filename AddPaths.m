function AddPaths()
%%%%%%%%%%%%%%%%%%%%%

disp('ADDING PATHS (2 subfolders to path)...')

addpath(genpath([pwd '\matlab']));
addpath(genpath([pwd '\mex']));
fprintf('\n');
