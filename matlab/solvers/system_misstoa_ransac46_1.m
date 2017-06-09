function [r,s,inliers,res,jac]=system_misstoa_ransac46_1(d,sys);
%
%
dobundle = 1;
if nargin < 2
    sys.ransac_threshold1 = 0.01;
    sys.ransac_k1 = 10;
end

% Overwrite parameters
if 1,
    sys.ransac_threshold1 = 0.01;
    sys.ransac_k1 = 10;
end


[r,s,inliers]=toa_3D_4n_ransac_orig(d,sys);
%sum(inliers(:))
if sum(inliers(:))>23,
    [r,s,res,jac]=toa_3D_bundle(d,r,s,inliers);
    [r,s]=toa_normalise(r,s);
end
%y=toa_trilateration(d,x);

%[m,n]=size(d); % Should really be nr of rec and send that we
% optimize over
%stdest = sqrt( (res'*res)/(length(res)-(3*m+3*n-6)) );




