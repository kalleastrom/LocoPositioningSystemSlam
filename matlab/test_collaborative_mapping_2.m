%% Collaborative mapping trial

addpath solvers

%% Generate a few measurement runs
% Output: runs
N = 2;

m = 6;
n = 200;
rdim = 2*m-3;
rsel = 1:(2*m);
rsel([1 2 4])=[];
sigma = 0.0001;
r = 10*rand(2,m);
s = 10*rand(2,n);
[r,s]=toa_normalise(r,s);
for k = 1:N,
    s = 10*rand(2,n); % Make new measurement positions for each run
    dtrue = toa_calc_d_from_xy(r,s);
    d = dtrue + sigma*randn(size(dtrue));
    runs(k).d = d;
    runs(k).rtrue = r;
    runs(k).strue = s;
end;

%% Run estimation on each run

for k = 1:N,
    % For each run bundle with ground truth as
    % initial estimate. This is a bit of cheating
    d = runs(k).d;
    [ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,runs(k).rtrue,runs(k).strue);
    % Normalize the coordinate system
    [ropt,sopt]=toa_normalise(ropt,sopt);
    dcalc = toa_calc_d_from_xy(ropt,sopt);
    res1 = (dcalc(:)-d(:));
    % Save the compressed parts (r,a,R)
    [m,n]=size(d);
    jaca = jacopt(:,1:(2*m));
    jacb = jacopt(:,(2*m+1):end);
    U = jaca'*jaca;
    W = jaca'*jacb;
    V = jacb'*jacb;
    dsdr = inv(V)*(W');
    jac_for_r = jaca+ jacb*dsdr;
    jac_for_r_normalized = jac_for_r(:,rsel);
    [qq,rr]=qr(full(jac_for_r_normalized));
    jac_qr = qq'*jac_for_r_normalized;
    jac_qr_small = jac_qr(1:rdim,1:rdim);
    % res_qr = qq'*resopt;
    % norm_res_qr = norm(res_qr);
    % Here comes the compressed result parts (a,R,r)
    a = norm(resopt);
    R = jac_qr_small;
    r = ropt;
    result(k).a = a;
    result(k).R = R;
    result(k).r = r;
    result(k).rvec = r(rsel)';
end;

%% Run estimation on all runs

dall = [];
sall = [];
for k = 1:N,
    % For each run bundle with ground truth as
    % initial estimate. This is a bit of cheating
    dall = [dall runs(k).d];
    sall = [sall runs(k).strue];
end
rall = runs(1).rtrue;

[ropt,sopt,resopt,jacopt]=toa_2D_bundle(dall,rall,sall);
% Normalize the coordinate system
[ropt,sopt]=toa_normalise(ropt,sopt);
dcalc = toa_calc_d_from_xy(ropt,sopt);
res1 = (dcalc(:)-dall(:));
% Save the compressed parts (r,a,R)
[m,n]=size(d);
jaca = jacopt(:,1:(2*m));
jacb = jacopt(:,(2*m+1):end);
U = jaca'*jaca;
W = jaca'*jacb;
V = jacb'*jacb;
dsdr = inv(V)*(W');
jac_for_r = jaca+ jacb*dsdr;
jac_for_r_normalized = jac_for_r(:,rsel);
[qq,rr]=qr(full(jac_for_r_normalized));
jac_qr = qq'*jac_for_r_normalized;
jac_qr_small = jac_qr(1:rdim,1:rdim);
% res_qr = qq'*resopt;
% norm_res_qr = norm(res_qr);
% Here comes the compressed result parts (a,R,r)
a = norm(resopt);
R = jac_qr_small;
r = ropt;
resultall.a = a;
resultall.R = R;
resultall.r = r;
resultall.rvec = r(rsel)';



%% What was the approximation again

% norm( res(rnew) ) \approx sqrt(a^2 + norm( R*(rnew-r))^2 )


%% Example 1. Estimate r from all compressed results

% Merge result 1 and 2
% minimize res_1^2 + res_2^2
% = a_1^2 + a_2^2 + R

M = [result(1).R;result(2).R];
b = [result(1).R*result(1).rvec;result(2).R*result(2).rvec];
r = M\b;

%% Also compress this result

w = M*r-b
aextra = w'*w
atot = sqrt(result(1).a^2 + result(2).a^2 + aextra);
[Qtot,Rtot]=qr(M);
Rtot = Rtot(1:rdim,1:rdim);

totresult.r = r;
totresult.R = Rtot;
totresult.a = atot;

%% Hyptesprövning, resdiuals should be 
%(m*n-parametrar)*sigma^2


% How well did this go?

[result(1).a^2 (m*n-2*n-rdim)*sigma^2]
[result(2).a^2 (m*n-2*n-rdim)*sigma^2]
[aextra (rdim)*sigma^2]
[resultall.a^2-result(1).a^2-result(2).a^2 aextra (rdim)*sigma^2]




