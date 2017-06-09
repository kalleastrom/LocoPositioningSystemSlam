%% Set up parameters for a set of experiments

%% Set parameters for experiment in ex
% make experimentlist
ex.NN = 10;  % Number of experiments
ex.m = 4;   % Number of microphones
ex.n = 25;   % Number of sounds
ex.outlier_ratio = 0.05; % Percentage of outlier measurements
% Each sound has m measurements, so this
% means that the chance that at least one of
% the measurements is wrong is 1-(1-outlier_ratio)^3
ex.missing_ratio = 0.00; % Percentage of missing measurements
ex.sigma = 0.001;    % Standard deviation of inlier errors
ex.outlier_low = 0.7; % Minimum errors for outliers
ex.outlier_high = 1.2;% Maximum errors for outliers
ex.outlier_n = round(ex.outlier_ratio*ex.m*ex.n);
ex.missing_n = round(ex.missing_ratio*ex.m*ex.n);
offset = 0;

clear exlist
clear exp_description

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increasing outlier ratio and noise simultaneously
llist = [0.00:0.02:0.10];
llist2 = llist/100;
ids = 1:length(llist);
for kk = ids;
    exlist(offset+kk)=ex;
    exlist(offset+kk).outlier_ratio = llist(kk);
    exlist(offset+kk).sigma = llist2(kk);
    exlist(offset+kk).outlier_n = round(exlist(offset+kk).outlier_ratio*ex.m*ex.n);
end;
exp_description(1).expid = ids + offset;
exp_description(1).variable = 'Outlier Ratio';
exp_description(1).variablefilename = 'Outlier_Ratio';
exp_description(1).variablevalue = llist;
offset = offset+length(llist);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increasing outlier ratio
% llist = [0.05:0.05:0.10];
% ids = 1:length(llist);
% for kk = ids;
%     exlist(offset+kk)=ex;
%     exlist(offset+kk).outlier_ratio = llist(kk);
% end;
% exp_description(1).expid = ids + offset;
% exp_description(1).variable = 'Outlier Ratio';
% exp_description(1).variablefilename = 'Outlier_Ratio';
% exp_description(1).variablevalue = llist;
% offset = offset+length(llist);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increasing missing data ratio
% llist = [0:0.04:0.24];
% ids = 1:length(llist);
% for kk = ids;
%     exlist(offset+kk)=ex;
%     exlist(offset+kk).missing_ratio = llist(kk);
% end;
% exp_description(2).expid = ids + offset;
% exp_description(2).variable = 'Missing Data Ratio';
% exp_description(2).variablefilename = 'Missing_Data_Ratio';
% exp_description(2).variablevalue = llist;
% offset = offset+length(llist);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % Increasing noise level
% llist = [0 0.001 0.0025 0.005 0.01 0.02 0.05];
% ids = 1:length(llist);
% for kk = ids;
%     exlist(offset+kk)=ex;
%     exlist(offset+kk).sigma = llist(kk);
% end;
% exp_description(3).expid = ids + offset;
% exp_description(3).variable = 'Inlier noise level';
% exp_description(3).variablefilename = 'Sigma';
% exp_description(3).variablevalue = llist;
% offset = offset+length(llist);
% 
% for kk = 1:length(exlist);
%     ex = exlist(kk);
%     exlist(kk).outlier_n = round(ex.outlier_ratio*ex.m*ex.n);
%     exlist(kk).missing_n = round(ex.missing_ratio*ex.m*ex.n);
% end
