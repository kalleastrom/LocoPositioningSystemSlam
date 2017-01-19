% fix a few path related things

thispath = pwd;
cd ..
rootpath = pwd;
cd(thispath);

addpath([rootpath filesep 'data']);
addpath([rootpath filesep 'matlab']);
addpath([rootpath filesep 'matlab' filesep 'solvers']);
addpath([rootpath filesep 'matlab' filesep 'simulations']);
addpath([rootpath filesep 'matlab' filesep 'tools']);

systemsettings.rootpath = rootpath;
