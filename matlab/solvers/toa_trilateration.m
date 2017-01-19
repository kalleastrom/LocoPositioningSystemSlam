function [y,inliers]=toa_trilateration(d,x,y0,index,inliers,ransac_k,ransac_tol);
% y=toa_trilateration_one_bundle(d,x,y0,index,inliers)
%



[m,n]=size(d);
[x_dim,xm]=size(x);

if nargin<3,
    y0 = NaN*ones(x_dim,n);
end

if nargin<4,
    index = find(ones(1,n));
end;

if nargin<5,
    inliers = isfinite(d);
end;

if nargin<6,
    ransac_k = 20;
end;

if nargin<7,
    ransac_tol = 0.2;
end;

y = y0;

for jj = index,
    %keyboard;
    [yy,inltmp,nr_inliers,err_rms] = toa_trilateration_one_ransac(d(:,jj),x,ransac_k,ransac_tol,find(inliers(:,jj)));
    inliers(:,jj)=inltmp(:);
    yy = real(yy);
    yy = toa_trilateration_one_bundle(d(:,jj),x,yy,inliers(:,jj));
    y(:,jj) = yy;
end;