function [r,s,inliers,sol]=system_misstoa_2D_ransac3n_bundle(d,sys);
%
%
if nargin<2,
    sys.ransac_k1 = 20;
    sys.ransac_tol = 1;
    sys.xdim = 2;
end;

% Get an initial estimate for three of the rows using ransac
% We hope to get as many inlier columns as possible
[sol,maxnrinl]=toa_misstoa_2D_3n_ransac(d,sys);
% This essentially runs [r0,s0,inliers]=toa_2D_3n_ransac_v2(d,sys);

% Doa little bundling
if sys.dobundle,
    [sol,res0,res,jac]=toa_misstoa_2D_bundle(sol,d,sys);
end

for blubb = 1:2,
    
    % Extend to more rows
    [sol]=toa_misstoa_2D_more_rows(sol,d,sys);
    if sys.dobundle,
        [sol,res0,res,jac]=toa_misstoa_2D_bundle(sol,d,sys);
    end
    
    % Extend to more cols
    [sol]=toa_misstoa_2D_more_cols(sol,d,sys);
    if sys.dobundle,
        [sol,res0,res,jac]=toa_misstoa_2D_bundle(sol,d,sys);
    end
    
end;

if 0,
    [rall,sall,res,jac]=toa_2D_bundle(d,r,s,ones(size(d)));
    dall=toa_calc_d_from_xy(rall,sall);
    
end

% Reshuffle the rows and cols
r = NaN*ones(2,size(d,1));
s = NaN*ones(2,size(d,2));
r(:,sol.rows)=sol.r;
s(:,sol.cols)=sol.s;
inliers = sol.inlmatrix;

