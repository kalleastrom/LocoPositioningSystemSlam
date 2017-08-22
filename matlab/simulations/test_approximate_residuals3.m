%% Test combining several solutions using approximate residuals
fix_paths;
clear;
%% Simulate
m = 4;         % Number of receivers
n = 50;       % Number of senders
sigma = 0.1; % Standard deviation measurement noise
n_batches = 10;

r = 10*rand(2,m);
s = 10*rand(2,n*n_batches);
% s = sortrows(s',1)'; % Sort to increase difference between batches.
[r,s]=toa_normalise(r,s);
dtrue = toa_calc_d_from_xy(r,s);
d = dtrue + sigma*randn(size(dtrue));
res0 = (dtrue(:)-d(:));

s_batch = reshape(s, 2, n, n_batches);
d_batch = reshape(d, m, n, n_batches);

%% Optimize all at once
[ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,r,s);

[ropt_all,sopt_all, fix_mask]=toa_normalise(ropt,sopt);
dcalc = toa_calc_d_from_xy(ropt,sopt);
res_all = (dcalc(:)-d(:));

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

res_tot = combine_compressed( valid_batch_res );
%% Evaluate
ropt_diff = zeros(1,n_batches);
xlabels = cell(1,n_batches);
for bi=1:n_batches
    ropt_diff(bi) = norm(ropt_all - batch_res(bi).ropt);
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