function [sol,maxnrinl]=toa_misstoa_2D_3n_ransac(d,sys);

sol.rows = [];
sol.cols = [];
sol.inlmatrix = zeros(size(d));
sol.d = d;
rows=1;
cols0=1;
iter=0;
while (size(d(rows,cols0),2) < 3)*(iter < 100)
<<<<<<< .mine
cols0 = [];
while size(cols0,2)<10,
    rows = randperm(size(d,1),3);
    cols0 = find(all(isfinite(d(rows,:))));
end
||||||| .r1060
rows = randperm(size(d,1));
rows = rows(1:3);
cols0 = find(all(isfinite(d(rows,:))));

=======
rows = randperm(size(d,1));
rows = rows(1:3);
cols0 = find(all(isfinite(d(rows,:))));
iter=iter+1;
end
>>>>>>> .r1088
[r0,s0,inliers]=toa_2D_3n_ransac_v2(d(rows,cols0),sys);
cols1 = find(sum(inliers));
cols = cols0(cols1);

sol.rows = rows;
sol.cols = cols;
sol.inlmatrix(rows,cols)=ones(length(rows),length(cols));
sol.r = r0;
sol.s = s0(:,cols1);

maxnrinl = sum(sum(sol.inlmatrix));
