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


%% Try some extra bundling

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


%% Visualize results

d = data.d;
dcalc = toa_calc_d_from_xy(rtmp,stmp);
resm = dcalc-d;
inl3 = (abs(resm)<T);

figure(4);
subplot(4,1,1);
imagesc(min(d,3));
subplot(4,1,2);
imagesc(min(dcalc,3));
subplot(4,1,3);
imagesc(dcalc-d);
subplot(4,1,4);
imagesc(inl3);

figure(5);
clf;
plot(d');
hold on;
plot(dcalc');

dtmp = d;
dtmp(find(~inl3))=NaN*ones(1,sum(sum(~inl3)));
figure(6);
clf;
plot(dtmp');
hold on;
plot(dcalc');


figure(7);
clf;
hist(resm(:),100);

figure(8);
plot(dcalc(:),resm(:),'.')


dcalcgt = toa_calc_d_from_xy(data.Gtr,data.GTs);

