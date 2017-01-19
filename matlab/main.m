% main.m

%% read a database of benchmark examples
[data]=read_from_lpsdb(systemsettings,1); %

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

allres = benchmark_evaluate_systems_on_experiments(systems,data.db);

%% Save in an appropriate place

%save allres allres

%%

nexp = length(data_db);
nsys = length(systems);
nrep = exlist(1).NN;
OK_matrix = zeros(nsys,nexp);
rerr_tensor = zeros(nrep,nsys,nexp);
t_matrix = zeros(nsys,nexp);

for kk = 1:size(rerr_tensor,3);
    for ii = 1:length(systems);
        ok_matrix(ii,kk)=allres{kk}(ii).ok;
        mrerr_matrix(ii,kk)=allres{kk}(ii).rerr;
        rerr_tensor(:,ii,kk)=allres{kk}(ii).rerrs;
        t_matrix(ii,kk)=allres{kk}(ii).t;
    end
end
%mean(ok');

%test_plots(rerr_tensor,systems)

%results(1)
%results(2)
%results(3)

%% Plots and/or visualizations of the results

%% Calculation of inlier/outlier ROC-plots

expi = 1;
nexp = length(allres);
oneres = allres{expi};
nsys = length(oneres);
aroc = zeros(nsys,nexp);
hrate = zeros(nsys,nexp);
ts = zeros(nsys,nexp);
rerrs = zeros(nsys,nexp);

for expi = 1:length(allres);
    oneres = allres{expi};
    for sysi = 1:length(oneres);
        oneres = allres{expi};
        onesysres = oneres(sysi);
        %yfacit = onesysres.est_inl(:);
        [m,n]=size(data_db{expi}(1).d);
        yfacit = zeros(m,n,length(data_db{expi}));
        for kk = 1:length(data_db{expi});
            yfacit(:,:,kk)=data_db{expi}(kk).Iinl;
        end
        yvalue = onesysres.est_res(:);
        okid = find(isfinite(yvalue));
        yvalue = abs(yvalue(okid));
        yfacit = yfacit(okid);
        
        [fpr,tpr,area_roc]=calcroc(yvalue,yfacit,1);
        aroc(sysi,expi)=area_roc;
        %hrate(sysi,expi)=allres{expi}(sysi).ok;
        hrate(sysi,expi)=mean(allres{expi}(sysi).rerrsopt < 0.1 );
        ts(sysi,expi)=allres{expi}(sysi).t;
        rerrs(sysi,expi)=allres{expi}(sysi).rerr;
        
    end
end

rerrs
hrate

%% Make plots

selplots = [3 6 5 4 2];
plotstyles = {'g-','b--','r--','g--','m--'};
for ploti = 1:length(exp_description);
    figure(ploti); clf;
    variablename = exp_description(ploti).variable;
    vfname = exp_description(ploti).variablefilename;
    xdata = exp_description(ploti).variablevalue;
    expid = exp_description(ploti).expid;
    ydata = hrate(selplots,expid);
    %     hold off;
    %     for blubb = 1:length(selplots);
    %       ph = plot(xdata',ydata(blubb,:)',plotstyles{blubb},'LineWidth',5);
    %       hold on;
    %     end;
    ph = plot(xdata',ydata','LineWidth',5);
    for blubb = 2:length(selplots),
        set(ph(blubb),'LineStyle','--');
    end
    %set(ph,'LineSpec','g-b--r--y--m--');
    title(['Success ratio as a function of ' variablename]);
    xlabel(variablename);
    ylabel('Probability of estimate close to the ground truth');
    ll = legend(systemtexts(selplots));
    set(ll,'FontSize',18);
    printFig(['fig_hitrate_vs_' vfname]);
end

%%
if 0,
    for ploti = 1:length(exp_description);
        figure(ploti);
        vfname = exp_description(ploti).variablefilename;
        printFig(['fig_hitrate_vs_' vfname]);
    end
end;
