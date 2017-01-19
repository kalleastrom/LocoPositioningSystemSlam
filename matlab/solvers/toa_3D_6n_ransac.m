function [x1,y1,inliers1]=toa_3D_6n_ransac(d);
% [x1,y1,inliers1]=toa_3D_6n_ransac(d);
% RANSAC based solver for overdetermined 
% case TOA node calibration
%  6 receivers and n senders, where n >= 5
% Measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%
% In: 
%    d - 6xn matrix with (n>=6) distances 
% Out: 
%    x1 - 3x6 matrix
%    y1 - 3xn matrix
%    inliers1 - logial matrix with inliers. 
%
% Idea 
%  * select random 4 columns
%  * Solve minimal case (6x4)
%  * Test to see how many of the remaining points can be reconstructed
%    within certain accuracy

[m,n]=size(d);
if m~=6,
    error('input d is assumed to have 5 rows');
end
if n<5,
    error('input d is assumed to have at least 6 columns');
end

load toa_3D_46_settings.mat
ransac_k1 = 15;
ransac_tol = 0.2;
xdim = 3;
bestx = NaN;
besty = NaN;
bestinliers = zeros(size(d));

%d=toa_calc_d_from_xy(xtmp,ytmp);

maxnrinliers = 0;
for kk = 1:ransac_k1
    pp = randperm(n);
    pp = pp(1:4);
    [sols,stats] = toa_3D_46(d(:,pp),settings);
    if size(sols,2)>0,
        nrsols = size(sols.x,2);
        if nrsols >= 1,
            for ii = 1:nrsols,
                xtmp = sols.x{ii};
                ytmp = sols.y{ii};
                [xtmp,ytmp]=toa_3D_bundle(d(:,pp),xtmp,ytmp);
                [xtmp,ytmp]=toa_normalise(xtmp,ytmp);
                nr_inliers = 0;
                
                inliers_d = zeros(m,n);
                y2 = NaN*ones(xdim,n);
                for jj = 1:n,
                    [y20,inliers0,nr_inliers0,err_rms]= ...
                        toa_trilateration_one_ransac(d(:,jj),xtmp,5,ransac_tol);
                    if nr_inliers0>=4,
                        nr_inliers = nr_inliers+1;
                        y2(:,jj)=y20;
                        inliers_d(:,jj)=inliers0;
                    else
                        inliers_d(:,jj)=zeros(m,1);
                    end;
                end
                nr_inliers = sum(sum(inliers_d));
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

