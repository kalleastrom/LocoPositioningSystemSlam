function [sol]=toa_misstoa_2D_more_rows(sol,d,sys);

dcut = d(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
rcut = sol.r;
scut = sol.s;

% Which rows (receivers) have not been estimated
rowsleft = setdiff(1:size(d,1),sol.rows);

for row = rowsleft,
    dcut = d(row,sol.cols)'; % dcut needs to be a column vector;
    xcut = scut;
    try
        [y,inliers,nr_inliers,err_rms]=toa_trilateration_one_ransac(dcut,xcut, ...
            sys.ransac_k,sys.ransac_tol);
        % Make a decision on adding this solution
        if nr_inliers >= sys.min_inliers2,
            sol.rows = [sol.rows row];
            sol.inlmatrix(row,sol.cols)=inliers';
            sol.r = [sol.r y];
        end
    catch
    end
end


