function [x1,y1,inliers1]=toa_2D_3n_ransac_v2(d,settings);
% [x1,y1,inliers1]=toa_2D_3n_ransac(d);
% RANSAC based solver for overdetermined
% case TOA node calibration
%  3 receivers and n senders, where n >= 4
% Measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%
% In:
%    d - 3xn matrix with (n>=7) distances
% Out:
%    x1 - 2x3 matrix
%    y1 - 2xn matrix
%    inliers1 - logial matrix with inliers.
%
% Idea
%  * select random 3 columns
%  * Solve minimal case (3x3)
%  * Test to see how many of the remaining points can be reconstructed
%    within certain accuracy

[m,n]=size(d);

if nargin<2,
    settings.ransac_k1 = 20;
    settings.ransac_threshold = 0.01;
    settings.xdim = 2;
end

ransac_k1 = settings.ransac_k1;
ransac_tol = settings.ransac_threshold;
%xdim = settings.xdim;

maxnrinliers = 0;
for kk = 1:ransac_k1
    pp = randperm(n);
    pp = pp(1:3);
    try
    [xsols,ysols] = toa_2D_33(d(:,pp));
    end
    %keyboard;
    if size(xsols,2)>0,
        nrsols = size(xsols,2);
        if nrsols >= 1,
            for ii = 1:nrsols,
                xtmp = xsols{ii};
                ytmp = ysols{ii};
                %[xtmp,ytmp]=toa_2D_bundle(d(:,pp),xtmp,ytmp);
                [xtmp,ytmp]=toa_normalise(xtmp,ytmp);
                
                %keyboard;
                r = xtmp;
                %                 (r(:,2:3)')\(d(2:3,:).^2 - ones(2,1)*(d(1,:).^2)) - (sum(r(:,2:3).^2)'*ones(1,size(d,2)))
                %
                %
                %                 -2*(r(:,2:3)')*s
                %                 (d(2:3,:).^2 - ones(2,1)*(d(1,:).^2)) - (sum(r(:,2:3).^2)'*ones(1,size(d,2)))
                %                                 -2*(r(:,2:3)')*s
                y2 = (-2*(r(:,2:3)'))\((d(2:3,:).^2 - ones(2,1)*(d(1,:).^2)) - (sum(r(:,2:3).^2)'*ones(1,size(d,2))));
                
                dtmp = toa_calc_d_from_xy(r,y2);
                
                inliers_tmp = all( abs(d-dtmp)<settings.ransac_tol);
                nr_inliers = sum(inliers_tmp);
                inliers_d = ones(3,1)*inliers_tmp;
                
                %nr_inliers
                %keyboard;
                
                if nr_inliers>maxnrinliers,
                    maxnrinliers = nr_inliers;
                    bestx = xtmp;
                    besty = y2;
                    bestinliers = inliers_d;
                end
            end
        end
    end
end
x1 = bestx;
y1 = besty;
inliers1 = bestinliers;
