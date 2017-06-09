function [xopt,yopt,res,jac]=toa_3D_bundle_with_smoother(d,x,y,inliers,opts);
% [xopt,yopt]=toa_3D_bundle(d,x,y,inliers);
% bundle adjustment, non-linear minimization of
%  sum_ij ( d(i,j)-sqrt(sum( (x(:,i)-y(:,j)).^2 )) )^2
% + smoothness prior
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

[xopt,yopt,res,jac]=bundletoa(D,I,J,x,y,0,opts);

%%

function [xopt,yopt,res,jac]=bundletoa(D,I,J,xt,yt,debug,opts);

if nargin<6,
    debug = 1;
end;

%debug = 1;
for kkk = 1:10;
    %kkk
    %keyboard

    [res,jac]=calcresandjac(D,I,J,xt,yt,opts);
    %dz = -(jac\res);
    dz = -(jac'*jac+0.1*eye(size(jac,2)))\(jac'*res);
    %dz = -(jac(:,3:end))\res;
    %dz = [zeros(2,1);dz];
    %     [u,s,v]=svd(full(jac),0);
    %     u = u(:,1:(end-6));
    %     s = s(1:(end-6),1:(end-6));
    %     v = v(:,1:(end-6));
    %     dz = -v\s*u'*res;
    %dz=-pinv(full(jac))*res; %THIS VERSON works for all sizes of jac, even when we ahve more unkowns than equations. When we have more equations tha unknowns (which is usually the case) it handels the extra DOF, that equates to singlar values=0, automativally
    [xtn,ytn]=updatexy(xt,yt,dz);
    [res2,jac2]=calcresandjac(D,I,J,xtn,ytn,opts);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    kkkk = 1;
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
                [res2,jac2]=calcresandjac(D,I,J,xtn,ytn,opts);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        kkkk
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


function [res,jac]=calcresandjac(D,I,J,x,y,opts)


[res1,jac1]=calcresandjac_toa(D,I,J,x,y);
if ~isempty(opts.cc)
    [res2,jac2]=calcresandjac_cc(opts.cc,x,y,opts.lambdacc);
else
    res2=[];
    jac2=[];
end

% if ~isempty(gps)
%     [res3,jac3]=calcresandjac_2D_gps(gps,x,y,opts.lambdagps);
% else
res3=[];
jac3=[];
% end
% res = [res1;res2];
% jac = [jac1;jac2];
res = [res1;res2;res3];
jac = [jac1;jac2;jac3];


function [res,jac]=calcresandjac_toa(D,I,J,x,y);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd-D;
II = (1:length(I))';
JJ1 = (I-1)*3+1;
JJ2 = (I-1)*3+2;
JJ3 = (I-1)*3+3;
JJ4 = (J-1)*3+1+3*m;
JJ5 = (J-1)*3+2+3*m;
JJ6 = (J-1)*3+3+3*m;

VV1 = idd.*Vt(:,1);
VV2 = idd.*Vt(:,2);
VV3 = idd.*Vt(:,3);
VV4 = -idd.*Vt(:,1);
VV5 = -idd.*Vt(:,2);
VV6 = -idd.*Vt(:,3);

jac = sparse([II;II;II;II;II;II],[JJ1;JJ2;JJ3;JJ4;JJ5;JJ6],[VV1;VV2;VV3;VV4;VV5;VV6],nn,3*m+3*n);

function [res,jac]=calcresandjac_cc(cc,x,y,lambda)

j1 = cc(:,1);
j2 = cc(:,2);
j3 = cc(:,3);

res1 = y(1,j1) + y(1,j3) - 2*y(1,j2);
res2 = y(2,j1) + y(2,j3) - 2*y(2,j2);
res3 = y(3,j1) + y(3,j3) - 2*y(3,j2);
nn = size(cc,1);
m = size(x,2);
n = size(y,2);
res = [res1';res2';res3'];
res = lambda *res;
% res = [res;n?gra residualer f?r linjeavvikelse];
IIx = (1:(nn))'; %How many residuals there are for first coord
IIy = ((nn+1):(2*nn))'; %And for the second coord
IIz = ((2*nn+1):(3*nn))'; %And for the second coord
JJ1x = (j1-1)*3+1+3*m; %All the indices of jac cols that has used first coord of j1
JJ1y = (j1-1)*3+2+3*m; %second coord of j1
JJ1z = (j1-1)*3+3+3*m; %third coord of j1
JJ2x = (j2-1)*3+1+3*m; %...
JJ2y = (j2-1)*3+2+3*m;
JJ2z = (j2-1)*3+3+3*m;
JJ3x = (j3-1)*3+1+3*m;
JJ3y = (j3-1)*3+2+3*m;
JJ3z = (j3-1)*3+3+3*m;

VV  = ones(size(JJ1x)); %derivativs are only ones, or minus 2*ones

jac = sparse([IIx;IIx;IIx;IIy;IIy;IIy;IIz;IIz;IIz], ...
    [JJ1x;JJ2x;JJ3x;JJ1y;JJ2y;JJ3y;JJ1z;JJ2z;JJ3z], ...
    [VV;-2*VV;VV;VV;-2*VV;VV;VV;-2*VV;VV],3*nn,3*m+3*n);
jac=jac*lambda;

function [xny,yny]=updatexy(x,y,dz);

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(3*m));
dz2 = dz((3*m+1):end);
xny = x + reshape(dz1,3,m);
yny = y + reshape(dz2,3,n);


