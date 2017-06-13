%% Collaborative mapping trial

addpath solvers

%% Generate a few measurement runs
% Output: runs
N = 5;

m = 4;
n = 200;
sigma = 0.001;
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
    [ropt,sopt,resopt,jacopt]=toa_2D_bundle(d,r,s);
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
    jac_for_r_normalized = jac_for_r(:,[3 5 6 7 8]);
    [qq,rr]=qr(full(jac_for_r_normalized));
    jac_qr = qq'*jac_for_r_normalized;
    jac_qr_small = jac_qr(1:5,1:5);
    % res_qr = qq'*resopt;
    % norm_res_qr = norm(res_qr);
    % Here comes the compressed result parts (a,R,r)
    a = norm(resopt);
    R = jac_qr_small;
    r = ropt;
    result(k).a = a;
    result(k).R = R;
    result(k).r = r;
    result(k).rvec = r([3 5 6 7 8])';
end;

%% What was the approximation again

% norm( res(rnew) ) \approx sqrt(a^2 + norm( R*(rnew-r))^2 )


%% Example 1. Estimate r from all compressed results

% Merge result 1 and 2
% minimize res_1^2 + res_2^2
% = a_1^2 + a_2^2 + R

M = [result(1).R;result(2).R];
b = [result(1).R*result(1).rvec;result(2).R*result(2).rvec];
r = M\b;
% Also compress this result

% How well did this go?


