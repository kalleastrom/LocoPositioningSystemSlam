function [xopt,yopt,res,jac]=toa_3D_bundle(d,x,y,inliers);
% [xopt,yopt]=toa_3D_bundle(d,x,y,inliers);
% bundle adjustment, non-linear minimization of
%  sum_ij ( d(i,j)-sqrt(sum( (x(:,i)-y(:,j)).^2 )) )^2
% Input: 
%   d - measurement matrix - size mxn
%   x - 3xm matrix initial estimate
%   y - 3xn matrix initial estimate
%   inliers - (optional) mxn logical matrix of inliers

if nargin<4,
    inliers = isfinite(d);
end
%keyboard;
[I,J,D]=find(inliers);
ind = sub2ind(size(d),I,J);
D = d(ind);

[xopt,yopt,res,jac]=bundletoa(D,I,J,x,y);
