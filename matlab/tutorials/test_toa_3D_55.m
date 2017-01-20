function test_toa_3D_55;

% test case 3D
options.int = 0;
options.normalise=1;
options.linear=0;
options.linearmotion=0;
options.planarmotion=0;
m = 5;
n = 5;
D = 3; % Dimensionality
options.twod = (D==2);

%addpath ..

[d,x,y] = simulate_toa(m,n,options);
[x,y]=toa_normalise(x,y);

load toa_3D_55_settings.mat
[sols,stats] = toa_3D_55(d,settings);

rmserr = zeros(1,size(sols.x,2));
for kk = 1:size(sols.x,2),
    xn = sols.x{kk};
    yn = sols.y{kk};
    [xn,yn]=toa_3D_bundle(d,xn,yn);
    [xn,yn]=toa_normalise(xn,yn);
    rmserr(kk) = sqrt( sum(sum( (xn-x).^2 )) + sum(sum( (yn-y).^2 )) );
end

[minrmserr,mini]=min(rmserr);
%log10(minrmserr)

xn = sols.x{mini};
yn = sols.y{mini};
[xn,yn]=toa_3D_bundle(d,xn,yn);
[xn,yn]=toa_normalise(xn,yn);

assertVectorsAlmostEqual(x,xn);
assertVectorsAlmostEqual(y,yn);
