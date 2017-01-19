function [sol,res0,res,jac]=toa_misstoa_2D_bundle(sol,d,sys);

dcut = d(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
rcut = sol.r;
scut = sol.s;

dcalc=toa_calc_d_from_xy(rcut,scut);
res0 = dcalc(find(icut))-dcut(find(icut));
[r1,s1,res,jac]=toa_2D_bundle(dcut,rcut,scut,icut);
sol.res = res;
sol.jac = jac;
