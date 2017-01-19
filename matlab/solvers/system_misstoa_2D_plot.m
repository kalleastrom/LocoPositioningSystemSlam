function dall = system_misstoa_2D_plot(d,sol);

%
rtmp = NaN*ones(2,size(d,1));
stmp = NaN*ones(2,size(d,2));
rtmp(:,sol.rows)=sol.r;
stmp(:,sol.cols)=sol.s;

if 0,
    [rtmp,stmp]=toa_2D_normalise(rtmp,stmp);
end;
dall=toa_calc_d_from_xy(rtmp,stmp);

%
figure(1);
% Plot of structure and motion
subplot(1,3,1);
hold off;
rita(rtmp,'*');
hold on
for kk = 1:size(rtmp,2);
    if isfinite(rtmp(1,kk));
        text(rtmp(1,kk)+0.1,rtmp(2,kk)+0.1,num2str(kk),'FontSize',20);
    end
end
rita(stmp,'-');
% Histogram of residuals
subplot(1,3,2);
hist(sol.res,50);
% image of residuals and outliers
subplot(1,3,3);
imagesc([d;dall.*sol.inlmatrix;(d-dall).*sol.inlmatrix]);

if 1,
    figure(2);
    plot(d');
    figure(3);
    plot(dall');
    figure(4);
    plot(d'-dall')
end


