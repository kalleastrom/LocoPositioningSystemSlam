function [xopt,yopt,res,jac]=toa_2D_bundle(d,x,y,inliers);
% [xopt,yopt]=toa_2D_bundle(d,x,y,inliers);
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

%%
%

function [xopt,yopt,res,jac]=bundletoa(D,I,J,xt,yt,debug);

if nargin<6,
    debug = 0;
end;

for kkk = 1:20;
    %kkk
    [res,jac]=calcresandjac(D,I,J,xt,yt);
    %dz = -(jac\res);
    dz = -(jac'*jac+0.1*speye(size(jac,2)))\(jac'*res);
    %     [u,s,v]=svd(full(jac),0);
    %     u = u(:,1:(end-6));
    %     s = s(1:(end-6),1:(end-6));
    %     v = v(:,1:(end-6));
    %     dz = -v\s*u'*res;
    %dz=-pinv(full(jac))*res; %THIS VERSON works for all sizes of jav, even when we ahve more unkowns than equations. When we have more equations tha unknowns (which is usually the case) it handels the extra DOF, that equates to singlar values=0, automativally
    [xtn,ytn]=updatexy(xt,yt,dz);
    [res2,jac2]=calcresandjac(D,I,J,xtn,ytn);
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
                [xtn,ytn]=updatexy(xt,yt,dz);
                [res2,jac2]=calcresandjac(D,I,J,xtn,ytn);
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
        keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        xt = xtn;
        yt = ytn;
    else
        if debug,
            disp([num2str(kkk) '  stalled']);
        end
    end
end;

xopt = xt;
yopt = yt;

function [res,jac]=calcresandjac(D,I,J,x,y);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd-D;
II = (1:length(I))';
JJ1 = (I-1)*2+1;
JJ2 = (I-1)*2+2;
JJ3 = (J-1)*2+1+2*m;
JJ4 = (J-1)*2+2+2*m;

VV1 = idd.*Vt(:,1);
VV2 = idd.*Vt(:,2);
VV3 = -idd.*Vt(:,1);
VV4 = -idd.*Vt(:,2);

jac = sparse([II;II;II;II],[JJ1;JJ2;JJ3;JJ4],[VV1;VV2;VV3;VV4],nn,2*m+2*n);


function [xny,yny]=updatexy(x,y,dz);

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(2*m));
dz2 = dz((2*m+1):end);
xny = x + reshape(dz1,2,m);
yny = y + reshape(dz2,2,n);



