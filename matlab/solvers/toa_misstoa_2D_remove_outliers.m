function [sol]=toa_misstoa_2D_remove_outliers(sol,d,sys,T);

dcut = d(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
rcut = sol.r;
scut = sol.s;

dcalc=toa_calc_d_from_xy(rcut,scut);
res0 = dcalc(find(icut))-dcut(find(icut));
badid = find( abs(res0)>T);

tmp = find(icut);
badid2 = tmp(badid);

icut(badid2)=zeros(length(badid2),1);
sol.inlmatrix(sol.rows,sol.cols)=icut;
