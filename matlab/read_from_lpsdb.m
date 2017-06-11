function [data]=read_from_lpsdb(systemsettings,id,full_dirname)
%[a,settings]=read_from_sfsdb(tdoasystemsettings,1);
if nargin < 3
    full_dirname = '';
end

switch id(1),
    case 1
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'lpsdb' filesep 'data_db_20170119_20170119.mat']);
        data.db = tmp.data_db;
        load exp_list_20150924B;
        data.exlist = exlist;
        data.exp_description = exp_description;
        data.anchorDim = 3; % The dimensionality of the span of the anchor positions
        data.bitcrazeDim = 3; % The dimensionality of the span of the bitcraze positions
    case 2
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'lpsdb' filesep 'data_db_20170119_20170119.mat']);
        data.db = tmp.data_db;
        for i = 1:length(data.db)
            data.db{i}=data.db{i}(1:3);
        end
        load exp_list_20150924B;
        data.exlist = exlist;
        data.exp_description = exp_description;
        data.db = data.db([1:7]);
        data.exp_description = data.exp_description(1);
        data.exlist = data.exlist(1:7);        
        data.anchorDim = 3; % The dimensionality of the span of the anchor positions
        data.bitcrazeDim = 3; % The dimensionality of the span of the bitcraze positions
    case 3
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data4.mat']);
        data = tmp.data;
    case 4
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data5.mat']);
        data = tmp.data;
    case 5
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data6.mat']);
        data = tmp.data;
    case 6
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data7.mat']);
        data = tmp.data;
    case 7
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data9.mat']);
        data = tmp.data;
    case 8
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data10.mat']);
        data = tmp.data;
    case 9
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data11.mat']);
        data = tmp.data;
    case 10
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data12.mat']);
        data = tmp.data;
    case 11
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data13.mat']);
        data = tmp.data;
    otherwise
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'data' filesep 'data4.mat']);
        data = tmp.data;
end
        
end

