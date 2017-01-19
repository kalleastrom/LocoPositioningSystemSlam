function [res,jac]=calcresandjac(mdata,param);

Bhatv = [param.R;sum(param.U(:,mdata.I22).*param.V(:,mdata.J22))'];
d2vc = mdata.bhat2d2*Bhatv;
res = d2vc-mdata.d2vm;
%res = Bhatv;

% calculate d(Bhatv)/dz first
% first the param.R part this is easy
II1 = param.indzr';
JJ1 = II1;
dBdz = ones(size(II1));
% Then the U'V part. This is a little tricky
II2 = ((length(param.R) + 1):(length(param.R)+length(mdata.I22)))';
dBdU1 = param.V(1,mdata.J22)';
dBdU2 = param.V(2,mdata.J22)';
dBdU3 = param.V(3,mdata.J22)';
dBdV1 = param.U(1,mdata.I22)';
dBdV2 = param.U(2,mdata.I22)';
dBdV3 = param.U(3,mdata.I22)';
JJ21 = (mdata.I22-1)*3+1+param.nzr;
JJ22 = (mdata.I22-1)*3+2+param.nzr;
JJ23 = (mdata.I22-1)*3+3+param.nzr;
JJ24 = (mdata.J22-1)*3+1+3*(param.mm-1)+param.nzr;
JJ25 = (mdata.J22-1)*3+2+3*(param.mm-1)+param.nzr;
JJ26 = (mdata.J22-1)*3+3+3*(param.mm-1)+param.nzr;

jac = sparse([II1;II2;II2;II2;II2;II2;II2],[JJ1;JJ21;JJ22;JJ23;JJ24;JJ25;JJ26],[dBdz;dBdU1;dBdU2;dBdU3;dBdV1;dBdV2;dBdV3],length(Bhatv),3*(param.mm-1)+3*(param.nn-1)+param.nzr);
% then remove columns corresponding to U(1:3,1:3)
jac(:,(param.nzr+1):(param.nzr+9))=[];

jac = mdata.bhat2d2*jac;
