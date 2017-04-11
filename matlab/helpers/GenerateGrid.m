function [bounds,steps] = GenerateGrid(sourceW,sourceH,delta,minTx,maxTx,minTy,maxTy,minR,maxR,minS,maxS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bounds.tx = [minTx,maxTx];
bounds.ty = [minTy,maxTy];
bounds.r = [minR,maxR]; 
bounds.s = [minS,maxS]; 

steps.tx =  delta*sourceW/sqrt(2);
steps.ty = delta*sourceH/sqrt(2);
steps.r = delta*sqrt(2);
steps.s = delta/sqrt(2);

