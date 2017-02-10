%% Tutorial for the TOA 6,4 problem
% Generate data matrix d
% New solver, 2017-02-09

options.int = 0;
options.normalise=1;
options.linear=0;
options.linearmotion=0;
options.planarmotion=0;
m = 6;
n = 4;
D = 3; % Dimensionality
options.twod = (D==2);

[d,x,y] = simulate_toa(m,n,options);
[xgt,ygt]=toa_normalise(x,y);

%% Run solver, old version
%load ../solvers/toa_3D_46_settings.mat
%[sols,stats] = toa_3D_46(d,settings);
%% Run solver, old version
[sols] = toa_3D_46_red(d);
%sols = solver_toa_46_red(data0);


rmserr = zeros(1,size(sols.x,2));
for kk = 1:size(sols.x,2),
    xn = sols.x{kk};
    yn = sols.y{kk};
    [xn,yn]=toa_3D_bundle(d,xn,yn);
    [xn,yn]=toa_normalise(xn,yn);
    rmserr(kk) = sqrt( sum(sum( (xn-x).^2 )) + sum(sum( (yn-y).^2 )) );
end

%% There are up to 38 solutions. Find the right one

[minrmserr,mini]=min(rmserr);

%% Plot the best of the 38 solutions vs the ground truth one

xn = sols.x{mini};
yn = sols.y{mini};
[xn,yn]=toa_3D_bundle(d,xn,yn);
[xn,yn]=toa_normalise(xn,yn);

rmserr = sqrt(sum(sum(([x y]-[xn yn]).^2 )))/30;

%%
figure(1); clf;
rita3(x,'g+');
hold on
rita3(y,'r+');
rita3(xn,'go');
rita3(yn,'ro');
legend({'Ground Truth Receivers','Ground Truth Senders','Estimated Receivers','Estimated Senders'});
title(['RMS error: ' num2str(rmserr)]);
