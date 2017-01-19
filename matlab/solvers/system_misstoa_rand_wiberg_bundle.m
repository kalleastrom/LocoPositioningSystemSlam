function [r,s,inliers]=system_misstoa_rand_wiberg_bundle(d)


rr = 3;
[n,m]=size(d);

inliers = isfinite(d);
d(~inliers)=0;
colsum = sum(inliers,1);
rowsum = sum(inliers,2);
[~,idc]=max(colsum);
[~,idr]=max(rowsum);

d(~isfinite(d(:,idc)),idc)=rand;
d(idr,~isfinite(d(idr,:)))=rand;

nidc = setdiff(1:m,idc);
nidr = setdiff(1:n,idr);

idrows = [idr nidr];
idcols = [idc nidc];

d2 = d(idrows,idcols).^2;
cl = compactionmatrix(n);
cr = compactionmatrix(m);

Y = cl*d2*cr';
W = inliers(nidr,nidc);

[~,~,Vini]=svd(Y);

[U,V] = damped_wiberg_new(Y,W,rr,Vini(:,1:3));



B = d2(idrows,idcols);
[cl,dl]=compactionmatrix(size(B,1));
[cr,dr]=compactionmatrix(size(B,2));
Bhat = dl*B*dr';
%Btilde = cl*B*cr';
%[u,s,v]=svd(Btilde);
%s(4:end,:)=zeros(size(s,1)-3,size(s,2));
%Btilde = u*s*v';

Btilde = U*V';
Bhat(2:end,2:end)=Btilde;
    


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
Hi = inv(H);
meig = min(eig(Hi));
if meig<=0,
    lille = 1e-10;
    Hi = Hi + (lille-meig)*eye(3);
end
L = chol(Hi);    
%L = chol(inv(H));
r00 = inv(L')*(xt);
s00 = L*(yt+repmat(b,1,size(Bhat,2)));

% Reorder the points

r0 = zeros(3,size(d,1));
s0 = zeros(3,size(d,2));

r0(:,idrows)=r00;
s0(:,idcols)=s00;

toa_calc_d_from_xy(r0,s0);

%inliers = sol.inlmatrix==1;

% Final bundle among the inliers

[r1,s1,res,jac]=toa_3D_bundle(d,r0,s0,inliers);
[r,s]=toa_normalise(r1,s1);

%[m,n]=size(d); % Should really be nr of rec and send that we
% optimize over
%stdest = sqrt( (res'*res)/(length(res)-(3*m+3*n-6)) );


