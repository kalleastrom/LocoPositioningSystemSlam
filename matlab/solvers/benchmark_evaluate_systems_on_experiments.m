function allres = benchmark_evaluate_systems_on_experiments(systems,data_db);
% allres = benchmark_evaluate_systems_on_experiments(systems,data_db);
%


%% Evaluate each system on each dataset in the benchmark

%systems = {'system_misstoa_rand_bundle','system_misstoa_ransac_bundle'};

clear rerr ok
for kk = 1:length(data_db);
    data = data_db{kk};
    NN = length(data);
    [kk length(data_db)]
    %ex = exlist{kk};
    for ii = 1:length(systems);
        ttmps = zeros(1,NN);
        rerrs = zeros(1,NN);
        oks = zeros(1,NN);
        [m,n]=size(data(1).d);
        est_inl = zeros(m,n,NN);
        est_res = zeros(m,n,NN);
        for jj = 1:NN;
            %disp([ii jj]);
            tic;
            %[rtmp,stmp,inltmp]=feval(systems{ii},data(jj).d,data(jj).r,data(jj).s);
            [rtmp,stmp,inltmp]=feval(systems{ii},data(jj).d);
            % 
            dc = toa_calc_d_from_xy(rtmp,stmp);
            %restmp = dc-data(jj).dtrue;
            restmp = dc-data(jj).d;
            ttmps(jj) = toc;
            %rerr(ii,jj) = sqrt(sum( (rtmp(:)-data(jj).r(:)).^2 ));
            rerrs(jj) = sqrt(sum( (rtmp(:)-data(jj).r(:)).^2 )/size(rtmp,2));
            %rerropt(ii,jj)=sqrt(sum( (rtmp(:)-data(jj).r(:)).^2 ));
            rerrsopt(jj) = sqrt(sum( (rtmp(:)-data(jj).ropt(:)).^2 )/size(rtmp,2));
            oks(jj)=rerrs(jj)<0.5;
            oksopt(jj)=rerrsopt(jj)<0.1;
            est_inl(:,:,jj)=inltmp;
            est_res(:,:,jj)=restmp;
        end
        %keyboard;
        results(ii).tmps = ttmps;
        results(ii).t = mean(ttmps);
        results(ii).rerrs = rerrs;
        results(ii).rerr = mean(rerrs);
        results(ii).rerrsopt = rerrsopt;
        results(ii).rerropt = mean(rerrsopt);
        results(ii).system = systems{ii};
        results(ii).oks = oks;
        results(ii).ok = mean(oks);
        results(ii).oksopt = oksopt;
        results(ii).okopt = mean(oksopt);
        results(ii).est_inl = est_inl;
        results(ii).est_res = est_res;
    end;
    allres{kk}=results;
end

%ok = rerr < 0.01;

