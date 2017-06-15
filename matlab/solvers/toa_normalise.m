function [x1,y1,fixated_mask]=toa_normalise(x0,y0);
% [x1,y1]=toa_normalise(x0,y0);
% change coordinate system so that %
% x1(:,1) = 0
% x1(:,2:(dim+1)) is upper triangular with positive diagonal
%
xdim = size(x0,1);
m = size(x0,2);
n = size(y0,2);

% translation
t = -x0(:,1);
x = x0 + repmat(t,1,m);
y = y0 + repmat(t,1,n);
% rotation
[q,r]=qr(x(:,2:(1+xdim)));
x = q'*x;
y = q'*y;
% mirroring
M = diag(sign(diag(x(:,2:(1+xdim)))));
x1 = M*x;
y1 = M*y;

% Return fixation mask for x0
fixated = ones(size(x0));
fixated(:,2:end) = 1 - triu(fixated(:,2:end));
fixated_mask = logical(fixated(:));