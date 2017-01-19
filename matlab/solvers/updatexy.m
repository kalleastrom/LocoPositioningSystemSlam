function [param_new]=updatexy(param,dz);
%[param_new]=updatexy(param,dz);

param_new = param;
param_new.R = param_new.R + dz(param.indzr);
param_new.U(param.indu) = param_new.U(param.indu) + dz(param.indzu);
param_new.V(param.indv) = param_new.V(param.indv) + dz(param.indzv);
