function [ cres ] = compact_res( ropt, sopt, resopt, jacopt, fixated_mask)

rs_dim = size(ropt,1);
m = size(ropt,2);
n = size(sopt,2);
assert(n >= m);

jaca = jacopt(:,1:(rs_dim *m));
jacb = jacopt(:,(rs_dim *m+1):end);

% U = jaca'*jaca;
W = jaca'*jacb;
V = jacb'*jacb;

% figure(2);
% spy([U W;W' V])

dsdr = V\(W');

jac_for_r = jaca+ jacb*dsdr;
jac_for_r_normalized = jac_for_r(:,~fixated_mask);

[qq,rr]=qr(full(jac_for_r_normalized));
res_qr = qq'*resopt;
jac_qr = rr;

cres.norm_res = norm(res_qr);
cres.jac_qr_small = jac_qr(1:size(jac_qr,2),:);
cres.dsdr = dsdr;
cres.ropt_small = ropt(~fixated_mask);
cres.fixated = fixated_mask;


