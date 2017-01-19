function [sol,res0,res,jac]=toa_misstoa_2D_bundle_with_smoother(sol,d,sys);

i2i = NaN*ones(1,size(d,2));
i2i(sol.cols)=1:length(sol.cols);
ccny = sys.cc;
ccny(:)=i2i(ccny(:));

cccut = ccny;
ok = find(all(isfinite(ccny)'));
cccut = cccut(ok,:);
sys.cc = cccut;

dcut = d(sol.rows,sol.cols);
icut = sol.inlmatrix(sol.rows,sol.cols);
rcut = sol.r;
scut = sol.s;

dcalc=toa_calc_d_from_xy(rcut,scut);
res0 = dcalc(find(icut))-dcut(find(icut));
[r1,s1,res,jac]=toa_2D_bundle_with_smoother(dcut,rcut,scut,icut,sys);
sol.res = res;
sol.jac = jac;
sol.r = r1;
sol.s = s1;
