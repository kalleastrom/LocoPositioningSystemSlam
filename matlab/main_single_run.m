% As for now there are a few datasets in
% LocoPositioningSystemSlam/data/data/
% When running
%   read_from_lpsdb
% with data_nr = 3 to 11
% one can try a few datasets with at least partial ground truth.
%
% In the long run we aim to put the datasets separately on
% vision.maths.lth.se
%

% Before running, there are a few paths that needs to be set
%
fix_paths
addpath ../data/lpsdb   % datasets
addpath solvers         % Algorithms used for toa-calibration
addpath tools           % Other useful algorithms
addpath tutorials       % A few tutorials here

%%
data_nr = 4;   % Use example data data_nr
system_nr = 3; % Use system system_nr

%% Current methods
% bundle - non-linear least squares optimization
% bundle_l1 - non-linear l1 optimization of residuals

% 1 - Initialize receivers and senders randomly, then bundle
% 2 - Use compaction and svd on whole dataset to get initialization, then
% bundle
% 3 - Use method from paper (x) to initialize then bundle
% 4 -


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

%% read a database of benchmark examples
[data]=read_from_lpsdb(systemsettings,data_nr); % Use example data data_nr

%% solve
[rtmp,stmp,inltmp,res,jac]=feval(systems{system_nr},data.d);

%% Extra bundling
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
%
figure(1); clf;
subplot(3,1,1);
plot(data.GTs_resamp');
subplot(3,1,2);
plot(data.sopt');
subplot(3,1,3);
plot(stmp');

figure(2);
subplot(4,1,1);
imagesc(min(d,3));
title('measured distances d_ij');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');
subplot(4,1,2);
imagesc(min(dcalc,3));
title('Fitted distances d_ij');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');
subplot(4,1,3);
imagesc(dcalc-d);
title('Residuals in distances d_ij');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');
subplot(4,1,4);
imagesc(inl3);
title('Estimated as inliers I_ij');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');

%
figure(3);
clf;
subplot(2,1,1)
plot(d');
hold on;
plot(dcalc');
title('Measured and fitted distances');
ylabel('Distance (m)');
xlabel('Measurement nr (j)');
subplot(2,1,2);
dtmp = d;
dtmp(find(~inl3))=NaN*ones(1,sum(sum(~inl3)));
plot(dtmp');
hold on;
plot(dcalc');
title('Measured and fitted distances (inliers)');
ylabel('Distance (m)');
xlabel('Measurement nr (j)');

%
figure(4);
clf;
hist(resm(:),100);
title('Histogram of residuals (both inliers and outliers)');
%
figure(5);
plot(dcalc(:),resm(:),'.')
title('Residuals as a function of calculated distances');
xlabel('Estimated distance (m)');
ylabel('Distance residuals (m)');

%% Calculate distances from gt

dcalcgt = toa_calc_d_from_xy(data.Gtr,data.GTs);
dcalcgt_resamp = toa_calc_d_from_xy(data.Gtr,data.GTs_resamp);
dcalc = toa_calc_d_from_xy(rtmp,stmp);
%hack
%ojoj = find(isnan(dcalcgt(1,:)));
%dcalcgt(:,ojoj) = repmat(dcalcgt(:,ojoj(1)-1),1,length(ojoj));

figure(11);
subplot(3,1,1);
imagesc(min(dcalc,3));
title('Estimated distances');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');
subplot(3,1,2);
imagesc(min(dcalcgt,3));
title('Ground truth distances');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');
subplot(3,1,3);
imagesc(min(dcalcgt_resamp,3));
title('Ground truth distances');
ylabel('Anchor nr (i)');
xlabel('Measurement nr (j)');


figure(12); clf;
plot(dcalc');
hold on
plot(dcalcgt_resamp')
plot(d')

figure(13);
dcalcgt = toa_calc_d_from_xy(data.Gtr,data.GTs_resamp);
subplot(2,1,1);
plot(dcalc(:),resm(:),'.')
subplot(2,1,2);
plot(dcalcgt(:),dcalc(:)-dcalcgt(:),'.');

figure(15);
plot(dcalcgt(:),d(:)-dcalcgt(:),'.');



