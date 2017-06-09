%%

addpath solvers


%% Test of bundle and approximation of residuals

m = 4; 
n = 200;
sigma = 0.001;
r = 10*rand(2,m);
s = 10*rand(2,n);
[r,s]=toa_normalise(r,s);
dtrue = toa_calc_d_from_xy(r,s);
d = dtrue + sigma*randn(size(dtrue));

%%

dcalc = toa_calc_d_from_xy(r,s);
res0 = (dcalc(:)-d(:));

[ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,r,s);

[ropt,sopt]=toa_normalise(ropt,sopt);
dcalc = toa_calc_d_from_xy(ropt,sopt);
res1 = (dcalc(:)-d(:));

%%

figure(1);
spy(jacopt)
jaca = jacopt(:,1:(2*m));
jacb = jacopt(:,(2*m+1):end);

U = jaca'*jaca;
W = jaca'*jacb;
V = jacb'*jacb;

figure(2);
spy([U W;W' V])

dsdr = inv(V)*(W');

%%

% Study the error function for 
% maps of type 
% r = ropt + z*rdir;

% Generate random direction rdir in R^5
rdir = rand(5,1);
rdir = rdir/norm(rdir);
% Reshape as 2xm map
rdirm = reshape([0;0;rdir(1);0;rdir(2:end)],2,m);

zlist = 0:0.1:1;
errs = zeros(length(zlist),3);
for k = 1:length(zlist);
    z = zlist(k);
    dr_5 = z*rdir;
    dr = z*rdirm(:); % delta r (as vector)
    drm = z*rdirm;   % delta r (as matrix)
    this_r = ropt + z*rdirm; % perturb anchors
    % perturb the map optimally for this r
    ds = dsdr*dr;
    dsm = reshape(ds,2,n);
    this_s = sopt + dsm;
    % calculate residual
    dcalc = toa_calc_d_from_xy(this_r,this_s);
    this_res = (dcalc(:)-d(:));
    app_res = resopt + jacopt*[dr;ds];
    jac_for_r = jaca+ jacb*dsdr;
    jac_for_r_normalized = jac_for_r(:,[3 5 6 7 8]);
    app_res2 = resopt + jac_for_r_normalized*dr_5;
    [qq,rr]=qr(full(jac_for_r_normalized));
    res_qr = qq'*resopt;
    jac_qr = qq'*jac_for_r_normalized;
    norm_res_qr = norm(res_qr);
    jac_qr_small = jac_qr(1:5,1:5);
    if 0,
        figure(3);
        clf
        hold off;
        plot(this_res,'g');
        hold on;
        plot(app_res,'ro');
        plot(app_res2,'b+');
        pause;
    end
    [norm(app_res-this_res)/norm(this_res) norm(app_res2-this_res)/norm(this_res)]
    %pause;
    errs(k,:)=[norm(this_res) norm(app_res) norm(app_res2)];
end


figure(4); clf;
plot(zlist,errs);
hold on
plot(zlist,sqrt(norm_res_qr^2 + zlist.^2*norm(jac_qr_small*rdir)^2),'o');

result.ropt = ropt;
result.ropt_small = ropt([3 5 6 7 8])';
result.norm_res_qr = norm_res_qr;
result.jac_qr_small = jac_qr_small;

%%

size(jacopt)
size(r)
size(s)
result



