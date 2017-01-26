%% test case 3D - Generate data matrix d

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
[x,y]=toa_normalise(x,y);

%% Run solver
load toa_3D_46_settings.mat
[sols,stats] = toa_3D_46(d,settings);

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

xn = sols.x{mini};
yn = sols.y{mini};
[xn,yn]=toa_3D_bundle(d,xn,yn);
[xn,yn]=toa_normalise(xn,yn);

figure(1);
subplot(1,2,1);
display_solution(d,x,y);
title('Ground Truth Solution');
subplot(1,2,1);
display_solution(d,xn,yn);
title('Solution from minimal 4,6 solver');
