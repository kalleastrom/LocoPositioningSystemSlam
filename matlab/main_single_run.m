% Dowbnload data4.mat and put it in
%

% main.m
fix_paths
addpath ../data/lpsdb
addpath solvers
addpath tools
addpath tutorials

%%
data_nr = 3;   % Use example data data_nr
system_nr = 3; % Use system system_nr


%% read a database of benchmark examples
[data]=read_from_lpsdb(systemsettings,data_nr); % Use example data data_nr

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
[rtmp,stmp,inltmp,res,jac]=feval(systems{system_nr},data.d);



%% Visualize results

T = 0.2;
d = data.d;
dcalc = toa_calc_d_from_xy(rtmp,stmp);
resm = dcalc-d;
inl3 = (abs(resm)<T);

figure(4);
subplot(4,1,1);
imagesc(d);
subplot(4,1,2);
imagesc(dcalc);
subplot(4,1,3);
imagesc(dcalc-d);
subplot(4,1,4);
imagesc(inl3);

figure(5);
clf;
plot(d');
hold on;
plot(dcalc');

figure(6);
clf;
hist(resm(:),100);

