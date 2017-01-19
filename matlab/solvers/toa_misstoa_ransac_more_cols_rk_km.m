function [sol]=toa_misstoa_ransac_more_cols(d,sol,sys);

[m,n]=size(d);
d2 = d.^2;
trycols = setdiff(1:n,sol.cols);
[cl,dl]=compactionmatrix(length(sol.rows));
[u,s,v]=svd(sol.Bhat(2:end,2:end));
u = u(:,1:3);
%keyboard;

for ii = trycols,
    d2n = d2(sol.rows,ii);
    maxnrinl = 0;
    for kk =1:sys.ransac_k2;
        okrows = find(isfinite(d2n));
        tmp = randperm(length(okrows));
        if length(tmp)>=4,
            tryrows1 = okrows(tmp(1:4));
            %trycols2 = sol.cols(trycols1);
            % Equations are d2n = Dhat(1,:)*inv(dr') + x*[1 1;0 V] = zz+x*ZZ;
            zz = inv(dl)*sol.Bhat(:,1);
            ZZ = [ones(length(sol.rows),1) [zeros(1,3);u]];
            ZZ0 = [[1;zeros(length(sol.rows)-1,1)] [zeros(1,3);u]];
            xx = inv(ZZ(tryrows1,:))*(d2n(tryrows1,1)-zz(tryrows1,1));
            %[d2n(okcols);zz(okcols)+xx*ZZ(:,okcols)]
            inlids = find( abs(d2n(okrows) - (zz(okrows)+ZZ(okrows,:)*xx) ) < sys.ransac_threshold2 );
            if length(inlids)>maxnrinl,
                maxnrinl = length(inlids);
                clear tmpsol
                tmpsol.rows = [sol.rows(tryrows1)];
                tmpsol.col =  [ii];
                tmpsol.Bhatn = ZZ0*xx;
                tmpsol.inlrows = sol.rows(okrows(inlids));
            end
        end;
    end
    if maxnrinl>sys.min_inliers2;
        sol.cols = [sol.cols tmpsol.col];
        sol.inlmatrix(tmpsol.inlrows,tmpsol.col) = ones(length(tmpsol.inlrows),1);
        sol.Bhat = [sol.Bhat tmpsol.Bhatn];
        sol.dl = compactionmatrix(length(sol.cols));
    end
end
