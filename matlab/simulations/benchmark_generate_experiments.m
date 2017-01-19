function [data,data_db]=benchmark_generate_experiments(exlist);
% [data,data_db]=benchmark_generate_experiments(exlist);

%% Generate problem dataset
clear data data_db;
for exp_k = 1:length(exlist);
    ex = exlist(exp_k);
    for k = 1:ex.NN;
        r = diag([10 10 3])*rand(3,ex.m); % Random placement in room 10x10x3 meters
        s = diag([10 10 3])*rand(3,ex.n); % Random placement in room 10x10x3 meters       
        [r,s]=toa_normalise(r,s);          % Change coordinate system for comparison
        d = toa_calc_d_from_xy(r,s);
        dtrue = d;
        d = d + randn(size(d))*ex.sigma;
        tmp0 = randperm(ex.m*ex.n);
        tmp1 = tmp0(1:ex.outlier_n);
        Iinl = ones(ex.m,ex.n);
        Iinl(tmp1)=zeros(1,length(ex.outlier_n));
        err = rand(1,ex.outlier_n)*(ex.outlier_high-ex.outlier_low)+ex.outlier_low;
        err = err.*sign(randn(size(err)));
        d(tmp1)=d(tmp1)+err; % We should make sure that d is positive
        Ioutl = 1-Iinl;
        tmp2 = tmp0((1:ex.missing_n)+ex.outlier_n);
        d(tmp2)=NaN*ones(ex.missing_n,1);
        Imiss = zeros(ex.m,ex.n);
        Imiss(tmp2)=ones(1,length(ex.missing_n));
        Iinl(tmp2)=zeros(1,length(ex.outlier_n));
        [ropt,sopt,res,jac]=toa_3D_bundle(d,r,s,Iinl);
        [ropt,sopt]=toa_normalise(ropt,sopt);          % Change coordinate system for comparison

        % Here we could have a test to see if the problem
        % is well defined by looking at the jacobean.
        
        data(k).r = r;
        data(k).s = s;
        data(k).ropt = ropt;
        data(k).sopt = sopt;
        data(k).d = d;
        data(k).dtrue = dtrue;
        data(k).d2true = dtrue.^2;
        data(k).Iinl = Iinl;
        data(k).Ioutl = Ioutl;
        data(k).Imiss = Imiss;
    end
    data_db{exp_k}=data;
end;

