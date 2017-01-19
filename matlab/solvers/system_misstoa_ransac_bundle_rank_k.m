function [r,s,inliers]=system_misstoa_ransac_bundle_rank_k(d)
%
%

dobundle = 1;

sys.ransac_threshold = 0.2;
sys.ransac_k = 70;
sys.ransac_threshold2 = 0.2;
sys.ransac_k2 = 20;
sys.min_inliers2 = 8;
sys.rk = 2;

if 0,
    sys.ransac_threshold = 5;
    sys.ransac_k = 70;
    sys.ransac_threshold2 = 2;
    sys.ransac_k2 = 20;
    sys.min_inliers2 = 4;
    sys.rk = 2;
end

for kk = 1:1;
    [sol,maxnrinl]=toa_misstoa_ransac_k_plus_2_rows(d,sys);
    if dobundle,
        [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,sys.rk);
    end
    
    [sol]=toa_misstoa_ransac_more_rows_rk_k(d,sol,sys);
    if dobundle,
        [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,sys.rk);
    end
    
    [sol]=toa_misstoa_ransac_more_cols_rk_k(d,sol,sys);
    if dobundle,
        [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,sys.rk);
    end
    
    [sol]=toa_misstoa_ransac_more_rows_rk_k(d,sol,sys);
    if dobundle,
        [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,sys.rk);
    end
    
    [sol]=toa_misstoa_ransac_more_cols_rk_k(d,sol,sys);
    if dobundle,
        [sol,res0,res,d2calc]=toa_misstoa_bundle_rank_k(sol,d,sys.rk);
    end
    
end;

rk = sys.rk;

% Go from Bhat via H,b to r,s
Bhat = sol.Bhat;
D = sys.rk;
[u,s,v]=svd(Bhat(2:end,2:end));
xr = u(:,1:D)';
yr = s(1:D,1:D)*v(:,1:D)';
xtp = [zeros(D,1) xr];
yt =[zeros(D,1) yr];
xt = xtp/(-2);
% Going to assume that there are many rows
Bhatcol1 = Bhat(:,1);

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
for i = 2:size(xt,2);
    eqs(i-1) = (-2*xt(:,i)'*bv + xt(:,i)'*Cv*xt(:,i)) - Bhatcol1(i);
end

eqs_linear = eqs;
[cfm_linear,mons_linear] = polynomials2matrix(eqs_linear);
cfm_linear = cfm_linear ./ repmat(sqrt(sum(cfm_linear.^2, 2)), 1, size(cfm_linear, 2));
cfm_linear0=cfm_linear;
AA = cfm_linear0(:,1:(end-1));
bb = cfm_linear0(:,end);
zz0 = -pinv(AA)*bb;
%         zz0 = -AA\bb;

H = evaluate(Cv,[zz0;0]);
b = evaluate(bv,[zz0;0]);
if min(eig(H))>0,
    L = chol(inv(H));
else
    % hack used if H is not positive definite.
    mins = min(eig(H));
    H= H + (-mins+0.1)*eye(3);
    L = chol(inv(H));
end
r00 = inv(L')*(xt);
s00 = L*(yt+repmat(b,1,size(Bhat,2)));

% Reorder the points

r0 = zeros(D,size(d,1));
s0 = zeros(D,size(d,2));

r0(:,sol.rows)=r00;
s0(:,sol.cols)=s00;

toa_calc_d_from_xy(r0,s0);

inliers = sol.inlmatrix==1;

% Final bundle among the inliers

[r1,s1,res,jac]=toa_2D_bundle(d,r0,s0,inliers);
[r,s]=toa_normalise(r1,s1);

%[m,n]=size(d); % Should really be nr of rec and send that we
% optimize over
%stdest = sqrt( (res'*res)/(length(res)-(3*m+3*n-6)) );




