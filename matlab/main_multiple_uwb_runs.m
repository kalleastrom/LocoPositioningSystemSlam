% Dowbnload data4.mat and put it in
%

% main.m
fix_paths
addpath ../data/lpsdb
addpath solvers
addpath tools
addpath tutorials

%%
data_nrs = 3:11;   % Use example data data_nr
system_nr = 3; % Use system system_nr

%% Choose one or several systems to test.
systems = {...
    'system_misstoa_rand_bundle',...
    'system_misstoa_hack_bundle',...
    'system_misstoa_ransac_bundle',...
    'system_misstoa_rand_wiberg_bundle',...
    'system_misstoa_rand_bundle_l1',...
    'system_misstoa_hack_bundle_l1'...
    };
systemtexts = {...
    'Rand init l_2 opt',...
    'SVD init l_2 opt',...
    'Ransac + l_2 opt',...
    'l_2 opt using Wiberg alg',...
    'rand init + l_1 opt',...
    'SVD init + l_1 opt'...
    };

%%
for data_nr = data_nrs,
    
    % read a database of benchmark examples
    [data]=read_from_lpsdb(systemsettings,data_nr); % Use example data data_nr
    % solve
    [rtmp,stmp,inltmp,res,jac]=feval(systems{system_nr},data.d);
    % Extra bundling
    T = 0.3;
    rin = rtmp;
    sin = stmp;
    d = data.d;
    dcalc = toa_calc_d_from_xy(rin,sin);
    resm = dcalc-d;
    inlin = (abs(resm)<T);
    mid = 2:(size(d,2)-1);
    sys.cc = [mid-1;mid;mid+1]';
    sys.lambdacc = 10;
    [rut,sut,res,jac]=toa_3D_bundle_with_smoother(d,rin,sin,inlin,sys);
    rtmp = rut;
    stmp = sut;
    % Visualize results
    d = data.d;
    dcalc = toa_calc_d_from_xy(rtmp,stmp);
    resm = dcalc-d;
    inl3 = (abs(resm)<T);
    %
    figure(1); clf;
    subplot(3,1,1);
    plot(data.GTs_resamp');
    subplot(3,1,2);
    plot(data.sopt');
    subplot(3,1,3);
    plot(stmp');
    pause;
end;