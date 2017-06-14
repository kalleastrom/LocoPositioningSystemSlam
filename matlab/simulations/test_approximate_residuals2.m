%% Test of bundle and approximation of residuals
clear;
fix_paths;
%% Simulate
m = 4;         % Number of receivers
n = 200;       % Number of senders
sigma = 0.001; % Standard deviation measurement noise

r = 10*rand(2,m);
s = 10*rand(2,n);
[r,s]=toa_normalise(r,s);
dtrue = toa_calc_d_from_xy(r,s);
d = dtrue + sigma*randn(size(dtrue));
res0 = (dtrue(:)-d(:));


%% Optimize
[ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,r,s);

[ropt,sopt]=toa_normalise(ropt,sopt);
dcalc = toa_calc_d_from_xy(ropt,sopt);
res1 = (dcalc(:)-d(:));

%% Compact residual
cres = compact_res(ropt, sopt, resopt, jacopt);

%% Study the error function for 
% maps of type 
% r = ropt + z*rdir;
n_sim = 1000;
zlist = 0:0.1:1;
errs = zeros(n_sim,length(zlist),3);
for si = 1:n_sim

    % Generate random direction rdir
    rdir = zeros(m*2,1);
    rdir_small = rand(sum(~cres.fixated),1);
    rdir_small = rdir_small/norm(rdir_small);
    rdir(~cres.fixated) = rdir_small;
    % Reshape as 2xm map
    rdirm = reshape(rdir,2,m);

    zlist = 0:0.1:1;    
    for k = 1:length(zlist);
        z = zlist(k);
        dr = z*rdirm(:); % delta r (as vector)
        this_r = ropt + z*rdirm; % perturb anchors
        % perturb the map optimally for this r
        ds = cres.dsdr*dr;
        dsm = reshape(ds,2,n);
        this_s = sopt + dsm;
        % calculate residual
        dcalc = toa_calc_d_from_xy(this_r,this_s);
        this_res = (dcalc(:)-d(:));
        app_res = resopt + jacopt*[dr;ds];
        errs(si,k,1:2)=[norm(this_res) norm(app_res)];
    end

    errs(si,:,3) = sqrt(cres.norm_res^2 + zlist.^2*norm(cres.jac_qr_small*rdir_small)^2);
    if 0
        figure(11); clf;
        plot(zlist,squeeze(errs(si,:,:)));
        legend({'true res','app res','compact res'}, 'Location', 'NorthWest')
        pause;
    end
end

max_diff = abs(1- (errs(:,end,1)./errs(:,end,3)));
figure(12);clf;
hist(max_diff)

% figure(13);clf;
% subplot(2,1,1)
% plot((errs(:,:,1)),'.')
% subplot(2,1,2)
% plot((errs(:,:,3)),'.')


