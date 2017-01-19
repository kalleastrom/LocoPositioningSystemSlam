%% Set parameters for experiment in ex
disp('setup experiments');
if 1,
    benchmark_setup_experiment_list
    save exp_list_20170119 exlist exp_description
else
    load exp_list_20170119
end

%% Generate problem dataset
disp('Generate data');
if 1,
    [data,data_db]=benchmark_generate_experiments(exlist);
    save data_db_20170119_20170119 data data_db
else
    load data_db_20170119_20170119
end
