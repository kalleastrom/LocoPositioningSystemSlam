function [fpr,tpr,area_roc]=calcroc(yvalue,yfacit,makeplot);
% []=calcroc(yvalue,yfacit);

if nargin<3,
    makeplot = 0;
end;

[yss,yi]=sort(yvalue);
yscramble = yfacit(yi);
tmp1 = cumsum(~yscramble/(sum(~yfacit)+eps));
tmp2 = cumsum(yscramble/(sum(yfacit)+eps));
% Perhaps one should generate the convex hull of the curve

fpr = [0; tmp1; 1];
tpr = [0; tmp2; 1];
if makeplot,
    plot(fpr,tpr);
    xlabel('1-specificity FP/N');
    ylabel('sensitivity TP/P');
end

% calculate area under roc curve
% Assume that both tpr and fpr are increasing
ilo = 1:(length(tpr)-1);
ihi = ilo+1;
area_roc = sum( (tpr(ihi)-tpr(ilo)).*(1-(fpr(ihi)+fpr(ilo))/2) );
