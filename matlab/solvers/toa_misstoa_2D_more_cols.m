function [sol]=toa_misstoa_2D_more_cols(sol,d,sys);

dcut = d(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
rcut = sol.r;
scut = sol.s;

% Which rows (receivers) have not been estimated
colsleft = setdiff(1:size(d,2),sol.cols);

for col = colsleft,
    dcut = d(sol.rows,col); % dcut needs to be a column vector;
    xcut = rcut;
    try
        [y,inliers,nr_inliers,err_rms]=toa_trilateration_one_ransac(dcut,xcut, ...
            sys.ransac_k,sys.ransac_tol);
        % Make a decision on adding this solution
        if nr_inliers >= sys.min_inliers2,
            sol.cols = [sol.cols col];
            sol.inlmatrix(sol.rows,col)=inliers;
            sol.s = [sol.s y];
        end
    catch
    end
end




