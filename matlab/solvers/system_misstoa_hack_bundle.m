function [r,s,inliers]=system_misstoa_hack_bundle(d,r,s)
%
%

[m,n]=size(d);

if nargin==3,
    r0 = r;
    s0 = s;
else
    [r0,s0]=system_misstoa_hack_initialization(d);
end
inliers = isfinite(d);
[r1,s1,res,jac]=toa_3D_bundle(d,r0,s0,inliers);

[r,s]=toa_normalise(r1,s1);

