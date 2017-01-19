function [x1,y1,inliers1]=toa_3D_4n_ransac(d);
% [x1,y1,inliers1]=toa_3D_4n_ransac(d);
% RANSAC based solver for overdetermined 
% case TOA node calibration
%  4 receivers and n senders, where n >= 7
% Measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%
% In: 
%    d - 4xn matrix with (n>=7) distances 
% Out: 
%    x1 - 3x4 matrix
%    y1 - 3xn matrix
%    inliers1 - logial matrix with inliers. 
%
% Idea 
%  * select random 6 columns
%  * Solve minimal case (4x6)
%  * Test to see how many of the remaining points can be reconstructed
%    within certain accuracy

[m,n]=size(d);

load toa_3D_46_settings.mat
ransac_k1 = 20;
ransac_tol = 1.18;
ransac_tol = 0.2;
xdim = 3;

maxnrinliers = 0;
for kk = 1:ransac_k1
    pp = randperm(n);
    pp = pp(1:6);
    try % skip iteration for nan values
    [sols,stats] = toa_3D_46(d(:,pp)',settings)
    catch
        continue
    end
    if size(sols,2)>0,
        nrsols = size(sols.x,2);
        if nrsols >= 1,
            for ii = 1:nrsols,
                xtmp = sols.y{ii};
                ytmp = sols.x{ii};
                [xtmp,ytmp]=toa_3D_bundle(d(:,pp),xtmp,ytmp);
                [xtmp,ytmp]=toa_normalise(xtmp,ytmp);
                nr_inliers = 0;
                
                y2 = NaN*ones(xdim,n);
                for jj = 1:n,
                    try
                    [y20,inliers0,nr_inliers0,err_rms]= ...
                        toa_trilateration_one_ransac(d(:,jj),xtmp,5,ransac_tol);
                    if nr_inliers0==4,
                        nr_inliers = nr_inliers+1;
                        y2(:,jj)=y20;
                        inliers_d(:,jj)=inliers0;
                    else
                        inliers_d(:,jj)=zeros(m,1);
                    end;
                    end
                end
                if nr_inliers>maxnrinliers,
                    maxnrinliers = nr_inliers
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

% figure
plot3(x1(1,:),x1(2,:),x1(3,:),'ob')
hold on
plot3(y1(1,:),y1(2,:),y1(3,:),'r.')

