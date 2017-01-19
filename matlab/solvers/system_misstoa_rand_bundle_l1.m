function [r,s,inliers]=system_misstoa_rand_bundle_l1(d,r,s);
%
%

[m,n]=size(d);
inliers = isfinite(d);
[I,J,D]=find(inliers);
ind = sub2ind(size(d),I,J);
D = d(ind);

if nargin==3,
    r0 = r;
    s0 = s;
else
    if 1,
        r0 = randn(3,m);
        s0 = randn(3,n);
    else
        % Make a hack init
        % First fill in missing value using means of rows and columns
        [r0,s0]=system_misstoa_hack_initialization(d);
    end
end

%[r1,s1,res,jac]=toa_2D_bundle(d,r0,s0,inliers);
%[r1,s1,res,jac]=toa_2D_bundle_robust(d,r0,s0,inliers);
lambda = 1;
maxIterOuter = 10;
maxIterInner = 10;
[r1,s1,res]=bundleL1IteratedToaMics3DSounds3D(D, I, J, r0, s0, lambda, maxIterOuter, maxIterInner);

[r,s]=toa_normalise(r1,s1);

