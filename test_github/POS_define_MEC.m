function [CBM] = POS_define_MEC(CBM, ED, ErrorFP, ErrorHP, ErrorMin)
% ...

% tag measured fluxes
vm      = [CBM.v(ED.ind)];

% add pairs of slack variables
CBM.e1      = sdpvar(size(vm,1),1);
CBM.m1      = sdpvar(size(vm,1),1);
CBM.e2      = sdpvar(size(vm,1),1);
CBM.m2      = sdpvar(size(vm,1),1);

% represent uncertainty
biasFP = max(  ErrorMin,abs(ED.wm*ErrorFP));  
biasHP = max(2*ErrorMin,abs(ED.wm*ErrorHP));  
e2max   = biasFP;    
m2max   = biasFP;
alpha   = -log(0.1)./(biasHP'-e2max');  
beta    = -log(0.1)./(biasHP'-m2max');  

% define MEC constraints
CBM.CB = [CBM.CB, CBM.v(ED.ind)==ED.wm+CBM.e1-CBM.m1+CBM.e2-CBM.m2];
CBM.CB = [CBM.CB, CBM.e2<=e2max];
CBM.CB = [CBM.CB, CBM.m2<=m2max];
CBM.CB = [CBM.CB, CBM.e1>=0];
CBM.CB = [CBM.CB, CBM.m1>=0];
CBM.CB = [CBM.CB, CBM.e2>=0];
CBM.CB = [CBM.CB, CBM.m2>=0];

% define objective
CBM.J  = alpha*CBM.e1+beta*CBM.m1;