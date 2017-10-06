%% Test combining several solutions using approximate residuals
fix_paths;
clear;
%% Simulate
m = 4;         % Number of receivers
n = 50;       % Number of senders per batch
sigma = 0.1; % Standard deviation measurement noise
n_batches = 10;

r = 10*rand(2,m);
[r,~]=toa_normalise(r,r);
s = zeros(2,n*n_batches);
dtrue = zeros(m, n*n_batches);
d = zeros(m, n*n_batches);
s_batch = zeros(2, n, n_batches);
d_batch = zeros(m, n, n_batches);

for nb=1:n_batches
    this_s = 10*rand(2,n);
    this_r = r;
    if nb==1
        % Move one receiver
        %this_r(:,end) = this_r(:,end) + 5*sigma;
    end
    
    this_dtrue = toa_calc_d_from_xy(this_r,this_s);
    this_d = this_dtrue + sigma*randn(size(this_dtrue));
    
    s_batch(:,:,nb) = this_s; 
    d_batch(:,:,nb) = this_d;
    
    b_start = 1 + (nb-1)*n;
    b_end = b_start + n - 1;
    dtrue(:,b_start:b_end) = this_dtrue;
    d(:,b_start:b_end) = this_d;
    s(:,b_start:b_end) = this_s;
end
res0 = (dtrue(:)-d(:));

%% Optimize all at once
[ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,r,s);

[ropt_all,sopt_all, fix_mask]=toa_normalise(ropt,sopt);
dcalc = toa_calc_d_from_xy(ropt,sopt);
res_all = (dcalc(:)-d(:));
global_res_norm = norm(res_all);

%% Optimize batch-wise
figure(1);clf;
hold on
batch_res = struct();
ropt_plot = zeros(size(r,1), size(r,2), n_batches);
for bi=1:n_batches
    this_d = d_batch(:,:,bi);
    this_s = s_batch(:,:,bi);
    [ropt,sopt,resopt,jacopt]=toa_2D_bundle(this_d,r,this_s);

    [ropt,sopt]=toa_normalise(ropt,sopt);
    dcalc = toa_calc_d_from_xy(ropt,sopt);
    res = (dcalc(:)-this_d(:));
    batch_res(bi).ropt = ropt;
    batch_res(bi).sopt = sopt;
    batch_res(bi).res = res;
    batch_res(bi).cres = compact_res(ropt, sopt, resopt, jacopt, fix_mask);
    
    plot(ropt(1,:), ropt(2,:),'x')
    plot(sopt(1,:), sopt(2,:), '.')
    ropt_plot(:,:,bi) = ropt;
end
figure(2);clf;
subplot(2,1,1)
plot(squeeze(ropt_plot(1,:,:))')
title('Reciever position')
ylabel('x')
subplot(2,1,2)
plot(squeeze(ropt_plot(2,:,:))')
ylabel('y')
xlabel('Batch')
%% Combine batches
[ valid_batch_res,  valid_batch] = validate_map_batches( batch_res, sigma, m*n );
disp('Valid batch:')
valid_batch
res_tot = combine_compressed( valid_batch_res );
%% Evaluate
ropt_diff = zeros(1,n_batches);
batch_residual = zeros(1,n_batches);
xlabels = cell(1,n_batches);
for bi=1:n_batches
    ropt_diff(bi) = norm(ropt_all - batch_res(bi).ropt);
    cres = batch_res(bi).cres;
    batch_residual(bi) = cres.norm_res;
    xlabels{bi} = sprintf('Batch %i',bi);
end
approx_diff = norm(ropt_all - res_tot.ropt);
xdata = 1:length(ropt_diff);
figure(3);clf;
hold on
plot(xdata,ropt_diff)
plot(xdata,repmat(approx_diff,1,n_batches));
ylim([0,max(ropt_diff)*1.1])
xlabel('Batch')
ylabel('Diff to optimal bundle')
legend('Batch ropt diff norm', 'Fused ropt diff norm');

rdim = 2*m-3;
expected_batch_res = (m*n-2*n-rdim)*sigma^2;
%expected_fused_res = (rdim)*sigma^2;
best_possible_fused_res = (n_batches*m*n-2*n*n_batches-rdim)*sigma^2;
figure(4);clf;
hold on 
plot(xdata,batch_residual)
plot(xdata,repmat(expected_batch_res,1,n_batches));
plot(xdata,repmat(res_tot.norm_res2,1,n_batches));
plot(xdata,repmat(best_possible_fused_res,1,n_batches));
xlabel('Batch')
ylabel('Norm of residuals')
legend('Batch residual norm', 'Expected batch res', 'Fused residual norm', 'Best possible fused norm');
