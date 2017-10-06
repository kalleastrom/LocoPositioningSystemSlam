function [ res ] = combine_compressed( batch_list )
%COMBINE_COMPRESSED Combines maps with compressed residuals

n_batches = length(batch_list);
n_param = length(batch_list(1).cres.ropt_small);
A = zeros(n_param*n_batches,n_param);
b = zeros(n_param*n_batches,1);
fix_mask = batch_list(1).cres.fixated;

tot_residual = 0
for bi=1:n_batches   
    r2 = bi*n_param;
    r1 = r2 - n_param + 1;
    A(r1:r2,:) = batch_list(bi).cres.jac_qr_small;
    b(r1:r2,:) = batch_list(bi).cres.jac_qr_small*batch_list(bi).cres.ropt_small;
    tot_residual = tot_residual + batch_list(bi).cres.norm_res^2; 
end

ropt_approx = zeros(size(batch_list(1).ropt));
rsol = A\b;
w = A*rsol-b;
extra_residual2 = w'*w;
tot_norm_residual = sqrt(tot_residual + extra_residual2);

ropt_approx(~fix_mask) = rsol;

res.ropt = ropt_approx;
res.norm_res = norm(A*rsol - b);
res.norm_res2 = tot_norm_residual;

%Combine residual

end

