function [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,rk);
%
%

%% Initializing
[mm,nn]=size(sol.Bhat);
[cl,dl]=compactionmatrix(mm);
[cr,dr]=compactionmatrix(nn);
d2 = d.^2;
% Indices for ok elements in d2 after reordering
d2shuffle =d2(sol.rows,sol.cols);
ishuffle = sol.inlmatrix(sol.rows,sol.cols);
[I,J]=find(ishuffle==1);
V = (1:length(I))';
d2vm = d2shuffle(sub2ind(size(d2shuffle),I,J));
% Indices for the elements in Bhat that we need to calculate
ineededinBhat = ishuffle==1;
ineededinBhat(:,1)=ones(mm,1);
ineededinBhat(1,:)=ones(1,nn);
[I2,J2]=find(ineededinBhat);
% First reorder so that the first column and first row are
% first
ord1 = find(J2==1);
ord2 = find( (I2==1) & (J2>1) );
ord3 = find( (I2>1) & (J2>1) );
I2 = I2([ord1;ord2;ord3]);
J2 = J2([ord1;ord2;ord3]);
nspecial = length(ord1)+length(ord2);
nrest = length(ord3);
I22 = I2( (nspecial+1):(nspecial+nrest) )-1;
J22 = J2( (nspecial+1):(nspecial+nrest) )-1;
V2 = (1:length(I2))';
ij2v = sparse(I2,J2,V2,mm,nn,length(I2));

%% Generate sparse mapping from bhat-elements to d2 elements
II = V;
JJ = ij2v(sub2ind(size(ij2v),I,J));
%
sel1 = find(I>1);
II = [II;V(sel1)];
JJ = [JJ;ij2v(sub2ind(size(ij2v),ones(size(sel1)),J(sel1)))];
%
sel2 = find(J>1);
II = [II;V(sel2)];
JJ = [JJ;ij2v(sub2ind(size(ij2v),I(sel2),ones(size(sel2))))];
%
sel3 = find( (J>1) & (I>1) );
II = [II;V(sel3)];
JJ = [JJ;ij2v(sub2ind(size(ij2v),ones(size(sel3)),ones(size(sel3))))];
%
bhat2d2 = sparse(II,JJ,ones(size(II)),length(V),length(V2));

%% Extract optimization parameters (U,V,Z)
Bhat = sol.Bhat;
[u,s,v]=svd(Bhat(2:end,2:end));
U = u(:,1:rk)';
V = s(1:rk,1:rk)*v(:,1:rk)';
R = Bhat(sub2ind(size(Bhat),I2(1:nspecial),J2(1:nspecial)));
%
Uchange = ones(size(U));
Uchange(1:rk,1:rk)=zeros(rk,rk);
[iu,ju]=find(Uchange);
indu = sub2ind(size(U),iu,ju);
%
Vchange = ones(size(V));
[iv,jv]=find(Vchange);
indv = sub2ind(size(V),iv,jv);

nzr = length(R);
nzu = sum(sum(Uchange));
nzv = sum(sum(Vchange));

indzr = 1:nzr;
indzu = (nzr+1):(nzr+nzu);
indzv = (nzr+nzu+1):(nzr+nzu+nzv);

param.U = U;
param.V = V;
param.R = R;
param.indu = indu;
param.indv = indv;
param.indzr = indzr;
param.indzu = indzu;
param.indzv = indzv;
param.nzr = nzr;
param.mm = mm;
param.nn = nn;
mdata.d2vm = d2vm;
mdata.bhat2d2 = bhat2d2;
mdata.I22 = I22;
mdata.J22 = J22;
mdata.rk = rk;

%% Test to calculate res and jac from param

%BB = zeros(mm,nn);
%BB(2:end,2:end)=param.U'*param.V;
%BB(sub2ind(size(Bhat),I2(1:nspecial),J2(1:nspecial)))=param.Z;
%BBv = [param.R;sum(U(:,I22).*V(:,J22))'];
%d2vc = bhat2d2*BBv;
%[d2vm d2vc]
%res = d2vm-d2vc;

%%

% param0 = param;
% dzz = zeros(172,1);
% dzz = randn(172,1);
% litet = 0.0001;
% param1 = updatexy(param0,dzz*litet);
% [res0,jac0]=calcresandjac(mdata,param0);
% [res1,jac1]=calcresandjac(mdata,param1);
% 
%[(res1-res0)/litet jac0*dzz]
%%


debug = 0;
%%%%%%%%%%%%%%%%%%%

[res0,jac0]=calcresandjac_rk_k(mdata,param);
%keyboard;
for kkk = 1:5;
    %kkk
    %keyboard;
    [res,jac]=calcresandjac_rk_k(mdata,param);
    dz = -(jac\res);
    %dz = -(jac'*jac+0.001*eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    %nrparam = size(jac,2);
    %dof = nrparam-6;
    %u = u(:,1:dof);
    %s = s(1:dof,1:dof);
    %v = v(:,1:dof);
    %dz = -v*inv(s)*u'*res;
    [param_new]=updatexy(param,dz);
    [res2,jac2]=calcresandjac_rk_k(mdata,param_new);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) & (norm(res)<norm(res2)),
                dz = dz/2;
                [param_new]=updatexy(param,dz);
                [res2,jac2]=calcresandjac_rk_k(mdata,param_new);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        %keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        param = param_new;
    else
        %disp([num2str(kkk) '  stalled']);
    end
end;

param_opt = param;

BB = zeros(mm,nn);
BB(2:end,2:end)=param_opt.U'*param_opt.V;
BB(sub2ind(size(Bhat),I2(1:nspecial),J2(1:nspecial)))=param_opt.R;

d2calc = zeros(size(d));
d2calc_shuffle = inv(dl)*BB*inv(dr');
d2calc(sol.rows,sol.cols) = d2calc_shuffle;

sol.Bhat = BB;

% function [res,jac]=calcresandjac(data,param);
% 
% Bhatv = [param.Z;sum(param.U(:,mdata.I22).*param.V(:,mdata.J22))'];
% d2vc = mdata.bhat2d2*Bhatv;
% res = mdata.d2vm-d2vc;
% 
% % calculate d(Bhatv)/dz first
% % first the param.Z part this is easy
% II1 = param.indzr';
% JJ1 = II1;
% dBdz = ones(size(II1));
% % Then the U'V part. This is a little tricky
% II2 = ((length(param.Z) + 1):(length(param.Z)+length(mdata.I22)))';
% dBdU1 = param.V(1,mdata.J22)';
% dBdU2 = param.V(2,mdata.J22)';
% dBdU3 = param.V(3,mdata.J22)';
% dBdV1 = param.U(1,mdata.I22)';
% dBdV2 = param.U(2,mdata.I22)';
% dBdV3 = param.U(3,mdata.I22)';
% JJ21 = (I22-1)*3+1+param.nzr;
% JJ22 = (I22-1)*3+2+param.nzr;
% JJ23 = (I22-1)*3+3+param.nzr;
% JJ24 = (J22-1)*3+1+3*(mm-1)+param.nzr;
% JJ25 = (J22-1)*3+2+3*(mm-1)+param.nzr;
% JJ26 = (J22-1)*3+3+3*(mm-1)+param.nzr;
% 
% jac = sparse([II1;II2;II2;II2;II2;II2;II2],[JJ1;JJ21;JJ22;JJ23;JJ24;JJ25;JJ26],[dBdz;dBdU1;dBdU2;dBdU3;dBdV1;dBdV2;dBdV3],length(Bhatv),3*(param.mm-1)+3*(param.nn-1)+param.nzr);
% 
% 
% function [param_new]=updatexy(param,dz);
% %[param_new]=updatexy(param,dz);
% 
% param_new = param;
% param_new.R = param.R + dz(param.indzr);
% param_new.U(param.indu) = param_new.U(param.indu) + dz(param.indzu);
% param_new.V(param.indv) = param_new.V(param.indv) + dz(param.indzv);
% 
% 
