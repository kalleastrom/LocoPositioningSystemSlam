function [sol]=toa_misstoa_ransac_more_rows_rk_k(d,sol,sys);


[m,n]=size(d);
d2 = d.^2;
tryrows = setdiff(1:m,sol.rows);
[cr,dr]=compactionmatrix(length(sol.cols));
[u,s,v]=svd(sol.Bhat(2:end,2:end));
v = v(:,1:sys.rk)';
rk = sys.rk;

for ii = tryrows,
    d2n = d2(ii,sol.cols);
    maxnrinl = 0;
    for kk =1:sys.ransac_k2;
        okcols = find(isfinite(d2n));
        tmp = randperm(length(okcols));
        if length(tmp)>=(sys.rk+1),
            
            trycols1 = okcols(tmp(1:(sys.rk+1)));
            %trycols2 = sol.cols(trycols1);
            % Equations are d2n = Dhat(1,:)*inv(dr') + x*[1 1;0 V] = zz+x*ZZ;
            zz = sol.Bhat(1,:)*inv(dr');
            ZZ = [ones(1,length(sol.cols));zeros(rk,1) v];
            ZZ0 = [1 zeros(1,length(sol.cols)-1);zeros(rk,1) v];
            xx = (d2n(1,trycols1)-zz(1,trycols1))*inv(ZZ(:,trycols1));
            %[d2n(okcols);zz(okcols)+xx*ZZ(:,okcols)]
            inlids = find( abs(d2n(okcols) - (zz(okcols)+xx*ZZ(:,okcols)) ) < sys.ransac_threshold2 );
            if length(inlids)>maxnrinl,
                maxnrinl = length(inlids);
                clear tmpsol
                tmpsol.row = ii;
                tmpsol.cols = [sol.cols(trycols1)];
                tmpsol.Bhatn = xx*ZZ0;
                tmpsol.inlcols = sol.cols(okcols(inlids));
            end
        end
    end
    if maxnrinl>sys.min_inliers2;
        sol.rows = [sol.rows tmpsol.row];
        sol.inlmatrix(tmpsol.row,tmpsol.inlcols) = ones(1,length(tmpsol.inlcols));
        sol.Bhat = [sol.Bhat;tmpsol.Bhatn];
        sol.dl = compactionmatrix(length(sol.rows));
    end
end