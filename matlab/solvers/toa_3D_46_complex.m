function [sols,stats] = toa_3D_46_complex(d,settings,xt,yt)
%[sols,stats] = toa_3D_46(d,settings,xt,yt)
% Solver for minimal case TOA node calibration
%  4 senders and 6 receivers
%  measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%  Goal of solver is to reconstruct x and y from d.
%  There are 38 solutions counted with multiplicity and
%  complex solutions.
% In: 
%    d - 6x4 matrix with distances 
%    settings - settings for solver use
%               load toa_3D_46_settings
%    xt - 
%    yt - 
% Out: 
%    sols  - cell array with all real and feasible solvers
%    stats - statistics from solver
if ~isfield(settings,'n_rmeqs'); settings.n_rmeqs = 10; end;
if ~isfield(settings,'threshold'); settings.threshold = 1e-13; end;

%%
template = settings.template;
th       = settings.threshold;
n_rmeqs  = settings.n_rmeqs;
ids      = settings.ids;
eids     = settings.eids;

m = 6;
n = 4;
D = 3;

d2 = d.^2;
cr = [-ones(1,n-1);eye(n-1)];
cl = [-ones(1,m-1);eye(m-1)]';

if nargin < 4
    %%
    s = cl*d2*cr;
    [uu,ss,vv]=svd(s);
    
    % ss(1,1)/ss(3,3)
    %
    % dss = diag(ss);
    %
    xr11 = (uu(:,1:D)*ss(1:D,1:D))';
    yr1 = vv(:,1:D)';
    %
    % rss = sqrt(ss(1:D,1:D));
    % xr22 = (uu(:,1:D)*rss)';
    % yr2 = rss*vv(:,1:D)';
    
    xr33 = (uu(:,1:D))';
    yr3 = ss(1:D,1:D)*vv(:,1:D)';
    
    % sc   = ss(1,1)/ss(3,3);
    % ssc  = sqrt(sc);
    % ssm  = diag([1 1 ssc]);
    % ssi  = diag([1 1 1/ssc]);
    % xr44 = (uu(:,1:D)*ssi)';
    % yr4 = ssm*ss(1:D,1:D)*vv(:,1:D)';
    
    if (ss(1,1)/ss(3,3) > 10)
        xr11 = xr33;
        yr1  = yr3;
    end
    
    xtp = [zeros(D,1) xr11];
    yt =[zeros(D,1) yr1];
    xt = xtp/(-2);
% else
%     xtp = [zeros(D,1) xt'];
%     yt =[zeros(D,1) yt];
%     xt = xtp/(-2);
end
% [ss(1) ss(3,3)]

%% construct linear constraints
nrofunknowns = (D*(D+1))/2 + D + 1;
E = eye(nrofunknowns);
for i = 1:nrofunknowns;
    xv(i) = multipol(1,E(i,:)');
end
one  = multipol(1,zeros(nrofunknowns,1));
zero = multipol(0,zeros(nrofunknowns,1));
if D==3,
    Cv = [xv(1) xv(2) xv(3); xv(2) xv(4) xv(5); xv(3) xv(5) xv(6)];
elseif D==2,
    Cv = [xv(1) xv(2); xv(2) xv(3)];
end
% Lv = [xv(1) 0 0 ; xv(2) xv(3) 0; xv(4) xv(5) xv(6)];
if D==3,
    bv = [xv(7) xv(8) xv(9)]';
elseif D==2,
    bv = [xv(4) xv(5)]';
end

%%


%% here comes the linear equations
adjC = adjv(Cv);
detC = detv(Cv);
for i = 2:size(xt,2);
    eqs(i-1) = (-2*xt(:,i)'*bv + xt(:,i)'*Cv*xt(:,i)) - (d(i,1)^2-d(1,1)^2);
end
eqsn(1) =  bv'*adjC*bv - detC*d(1,1)^2;
for j = 2:size(yt,2);
    eqsn(j) = ( 2*bv'*adjC*yt(:,j) + yt(:,j)'*adjC*yt(:,j)) - detC*(d(1,j)^2-d(1,1)^2);
end
%         evaluate(eqs,gtall);
%         evaluate(eqsn,gtall);

%% eliminate some variables (xx(1:5))

eqs_linear = eqs;
[cfm_linear,mons_linear] = polynomials2matrix(eqs_linear);
cfm_linear = cfm_linear ./ repmat(sqrt(sum(cfm_linear.^2, 2)), 1, size(cfm_linear, 2));
cfm_linear0=cfm_linear;
AA = cfm_linear0(:,1:(end-1));
bb = cfm_linear0(:,end);
zz0 = -pinv(AA)*bb;
%         zz0 = -AA\bb;
[u,s,v]=svd(AA);
NN = v(:,(m):(nrofunknowns-1));
%         NN = v(:,[5,7:9]);


%% construct multipol equations
nrofunknowns2 = nrofunknowns-(m-1);
E = eye(nrofunknowns2);
for i = 1:(nrofunknowns2);
    zzv(i,1) = multipol(1,E(i,:)');
end
one  = multipol(1,zeros(nrofunknowns2,1));
zero = multipol(0,zeros(nrofunknowns2,1));

if nrofunknowns2>1,
    xv = (one*NN)*zzv(1:(nrofunknowns2-1))+one*zz0;
else
    xv = (one*zz0);
end
if D==3,
    Cv = [xv(1) xv(2) xv(3); xv(2) xv(4) xv(5); xv(3) xv(5) xv(6)];
elseif D==2,
    Cv = [xv(1) xv(2); xv(2) xv(3)];
end
% Lv = [xv(1) 0 0 ; xv(2) xv(3) 0; xv(4) xv(5) xv(6)];
if D==3,
    bv = [xv(7) xv(8) xv(9)]';
elseif D==2,
    bv = [xv(4) xv(5)]';
end


%% here comes the equations
% d(1,1) = bv'*inv(Cv)*bv
clear eqs;
adjC = adjv(Cv);
detC = detv(Cv);
eqs(1) =  bv'*adjC*bv - detC*d(1,1)^2;
eqsx10(1) =  bv'*adjC*bv - zzv(end)*d(1,1)^2;
for j = 2:size(yt,2);
    eqs(j) = ( 2*bv'*adjC*yt(:,j) + yt(:,j)'*adjC*yt(:,j)) - detC*(d(1,j)^2-d(1,1)^2);
    eqsx10(j) = ( 2*bv'*adjC*yt(:,j) + yt(:,j)'*adjC*yt(:,j)) - zzv(end)*(d(1,j)^2-d(1,1)^2);
end
eqs(n+1)=detC-zzv(end);
eqstry = [eqsx10 eqs(end)];
%% multiply with multmons
coeffs_eqs = coeffsof_plist2vec(eqstry);
cfm1   = zeros(size(template.T));
cfm1(ids) = coeffs_eqs(eids);
cfm1   = full(cfm1);
mons1  = template.mon;

mons1v = monvec2matrix(mons1);
id1 = find(mons1v(5,:)==0);
id2 = find(mons1v(5,:)>0);
cfm2 = cfm1(:,[id1 id2]);
mons2v = mons1v(:,[id1 id2]);

id1p = 1:length(id1);
id2p = (length(id1)+1):(length(id1)+length(id2));
[qq,rr,ee]=qr(cfm2(:,id1p));

drr = diag(rr);
drr = abs((drr/drr(1))) < th;
drr = drr(1:end-1) - drr(2:end);
kk  = find(drr < 0);

ff   = qq(:,(kk+1):end)'*cfm2(:,id2p);
gg   = repmat((sqrt(sum(ff.^2,2))),1,size(ff,2));
ff   = ff./gg;

% if ~exist('mon','var')
%     mons3v = mons2v-repmat([0 0 0 0 1]',1,size(mons2v,2));
%     eqs4 = cm2eqs(ff,mons3v(1:4,id2p));
% end

eqs00 = eqs(1:4)';
for kk = 1:4;
    mm = monomials(eqs(kk));
    cc = coeffs(eqs(kk));
    eqs00(kk) = multipol(cc,mm(1:4,:));
end

etest2  = eqs00;

settings.basis_size = 120;
settings.dim = 38;
settings.action_variable = 2;
settings.C = ff(n_rmeqs:end,:);
settings.mon = settings.mon2;
%%
[sols_short stats] = polysolve_C_mon(etest2, settings);

% find real solutions
sols_real = find_realsol(sols_short);
n_real    = size(sols_real,2);


%% recover C and b
sols_full =  NN*sols_real + repmat(zz0,1,n_real);
err =  inf;
kk = 1;
ok = 0;
for i = 1:size(sols_full,2);
    c = (sols_full(1:6,i));
    CC{i} = [c(1) c(2) c(3) ; c(2) c(4) c(5) ; c(3) c(5) c(6)];
    try
        l    = chol(inv(CC{i}));
        ok   = 1;
    catch % if not positive-definite, throw away
        %         disp(' C is not positive definite')
        l    = 0;
        continue;
    end
    C{kk} = CC{i};
    b{kk} = sols_full(end-2:end,i);
    L{kk} = l;
    solsvec{kk} = sols_full(:,i);
    kk = kk+1;
end

if ok
    %% reconstruct x, y
    err =  inf;
    % clear err;
    for kk = 1:length(L);
        x{kk} = inv(L{kk}')*(xt);
        y{kk} = L{kk}*(yt+repmat(b{kk},1,n));
        dd = sqrt( sum( (kron(ones(1,n),x{kk}) - kron(y{kk},ones(1,m)) ).^2 , 1 ) );
        dd = reshape(dd,m,n);
        err(kk) = norm(d - dd);
    end
    
    % choose the one with smallest residue
    bestid = find(err == min(err));
    xbest   = x{bestid};
    ybest   = y{bestid};
    solsvec_best  = solsvec{bestid};
    
    errbest = min(err);
    sols.C = C;
    sols.L = L;
    sols.b = b;
    sols.x = x;
    sols.y = y;
    sols.bestid = bestid;
    sols.n_real = n_real;
    stats.n_real = n_real;
    stats.n_good = length(L);
else
    sols = [];
    stats.n_good = 0;
    %     stats = [];
end

