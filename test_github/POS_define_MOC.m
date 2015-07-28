function [CBmodel] = POS_define_MOC(N, irrev)
% ...

% define fluxes
CBmodel.v = sdpvar(size(N,2),1);

% define constraints
CBmodel.CB = [N*CBmodel.v==0, diag([irrev])*CBmodel.v>=0];

% save data
CBmodel.Data.N = N;
CBmodel.Data.irrev = irrev;