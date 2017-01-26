function [data]=read_from_lpsdb(systemsettings,id,full_dirname)
%[a,settings]=read_from_sfsdb(tdoasystemsettings,1);
if nargin < 3
    full_dirname = '';
end

switch id(1),
    case 1
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'lpsdb' filesep 'data_db_20170119_20170119.mat']);
        data.db = tmp.data_db;
        data.anchorDim = 3; % The dimensionality of the span of the anchor positions
        data.bitcrazeDim = 3; % The dimensionality of the span of the bitcraze positions
    case 2
        tmp = load([systemsettings.rootpath filesep 'data' filesep 'lpsdb' filesep 'data_db_20170119_20170119.mat']);
        data.db = tmp.data_db;
        data.anchorDim = 3; % The dimensionality of the span of the anchor positions
        data.bitcrazeDim = 3; % The dimensionality of the span of the bitcraze positions
    otherwise
%         dirname = ensure_dirname(full_dirname,'sound files');
%         fileExtension = '.aiff';
%         sound_files = dir([dirname '*' fileExtension]);
%         settings.mm = length(sound_files);
%         settings.channels = 1:settings.mm;
%         
%         for i = 1:1
%             sound_file = sound_files(i);
%             ainfo = audioinfo(sound_file.name);
%             settings.sr = ainfo.SampleRate;
%             a = zeros(settings.mm,ainfo.TotalSamples);
%             a(i,:) = audioread(sound_file.name);
%         end
%         
%         for i = 2:settings.mm
%             sound_file = sound_files(i);
%             a(i,:) = audioread(sound_file.name);
%         end
end
        
end

