disp('COMPILING MEX FILES...')

cd mex

disp('==> compiling ''Configs2Affine_mex.cpp'' (1 out of 3)');
mex Configs2Affine_mex.cpp

disp('==> compiling ''CreateList_mex.cpp'' (2 out of 3)');
mex CreateList_mex.cpp

disp('==> compiling ''EvaluateConfigs_mex.cpp'' (3 out of 3)');
mex EvaluateConfigs_mex.cpp

disp('==> compiling ''EvaluateConfigsVectorized_mex.cpp'' (3 out of 3)');
mex EvaluateConfigsVectorized_mex.cpp

disp('==> DONE!');
fprintf('\n');
cd ..