function [r,s,inliers,res,jac]=system_misstoa_ransac46_4(d,sys);
%
%
%keyboard;
dobundle = 1;
if nargin < 2
    sys.ransac_threshold1 = 0.01;
    sys.ransac_k1 = 40;
    sys.solver = 'toa_3D_46_red';
end

% Overwrite parameters
if 1,
    sys.ransac_threshold1 = 0.01;
    sys.ransac_k1 = 40;
    sys.solver = 'toa_3D_46_red_v3';
end

%[sols] = toa_3D_64_red(d)
%solver46 = 'toa_3D_46_red';
[r,s,inliers]=toa_3D_4n_ransac_special2(d,sys);
%sum(inliers(:))
if sum(inliers(:))>23,
    [r,s,res,jac]=toa_3D_bundle(d,r,s,inliers);
    [r,s]=toa_normalise(r,s);
end
%y=toa_trilateration(d,x);

%[m,n]=size(d); % Should really be nr of rec and send that we
% optimize over
%stdest = sqrt( (res'*res)/(length(res)-(3*m+3*n-6)) );




