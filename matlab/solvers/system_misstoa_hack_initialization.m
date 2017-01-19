function [r,s,inliers]=system_misstoa_hack_initialization(d);

[m,n]=size(d);
inliers = isfinite(d);
dhack = d;
missing = find(~isfinite(dhack));
dhack(missing) = zeros(length(missing),1);
inliers = isfinite(d);
meanr = sum(dhack.*inliers)./sum(inliers);
meanc = sum((dhack.*inliers)')./sum(inliers');
fillin = (repmat(meanr,m,1)+repmat(meanc',1,n))/2;
dhack(missing)=fillin(missing);
% Form Bhat
[cl,dl]=compactionmatrix(m);
[cr,dr]=compactionmatrix(n);
Bhat = dl*(dhack.^2)*dr';

D = 3;
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

r0=r00;
s0=s00;
toa_calc_d_from_xy(r0,s0);

% Final bundle among the inliers

%[r1,s1,res,jac]=toa_3D_bundle(d,r0,s0,inliers);
[r,s]=toa_normalise(r0,s0);

%[m,n]=size(d); % Should really be nr of rec and send that we
% optimize over
%stdest = sqrt( (res'*res)/(length(res)-(3*m+3*n-6)) );




