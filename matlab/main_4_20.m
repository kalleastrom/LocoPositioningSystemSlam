% main.m

addpath ../data/lpsdb

%% Supress rank warnings

if 0,
    w = warning('query','last')
    id = w.identifier;
    warning('off',id);
end


%% read a database of benchmark examples
%[data]=read_from_lpsdb(systemsettings,3); %
%load data_db_4_25
if 1,
    benchmark_setup_experiment_list_4_20;
    [data,data_db]=benchmark_generate_experiments(exlist);
    % save data_db_4_25_short data_db exlist exp_description
else
    load data_db_4_25_short
end
%% Choose one or several systems to test.
systems = {...
    'system_misstoa_ransac46_1',...
    'system_misstoa_ransac46_2',...
    'system_misstoa_ransac46_3',...
    'system_misstoa_ransac46_4',...
    };
systemtexts = {...
    'Kuang et al ',...
    'New solver',...
    'New solver more iterations',...
    'New solver2 more iterations',...
    };

%%

% Make routines for benchmarking with or without smoothing
allres = benchmark_evaluate_systems_on_experiments(systems,data_db);


%%
% ii = 1;
% [rtmp,stmp,inltmp]=feval(systems{ii},data.d);


%% Save in an appropriate place

%save allres allres

%%

nexp = length(data_db);
nsys = length(systems);
%nrep = exlist(1).NN;
nrep = length(allres{1}(1).oks);
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

%%
variablename='Outlier percentage';

%% Make plots

selplots = [1 2 3 4];
plotstyles = {'r-','b--','y:','g:'};
plotstyles2 = {'-','--',':','-.'};
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
    ph = plot(100*xdata',100*ydata','LineWidth',5);
    for blubb = 2:length(selplots),
        set(ph(blubb),'LineStyle','--');
    end
    for blubb = 2:length(selplots),
        set(ph(blubb),'LineStyle',plotstyles2{blubb});
    end
    %set(ph,'LineSpec','g-b--r--y--m--');
    title(['Success ratio as a function of ' variablename],'FontSize',20);
    xlabel(variablename,'FontSize',20);
    ylabel('Success percentage','FontSize',20);
    ll = legend(systemtexts(selplots));
    set(ll,'FontSize',18);
    axis([0 10 0 100]);
%    printFig(['fig_hitrate_vs_' vfname]);
end


%%
if 0,
    for ploti = 1:length(exp_description);
        figure(ploti);
        vfname = exp_description(ploti).variablefilename;
        printFig(['fig_hitrate_vs_' vfname]);
    end
end;
