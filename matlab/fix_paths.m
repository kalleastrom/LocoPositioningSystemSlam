% fix a few path related things

% Restore path to avoid name clashes with other projects
restoredefaultpath;

thispath = pwd;
cd ..
rootpath = pwd;
cd ..
workspace = pwd;
cd(thispath);

addpath([rootpath filesep 'data']);
addpath([rootpath filesep 'matlab']);
addpath([rootpath filesep 'matlab' filesep 'solvers']);
addpath([rootpath filesep 'matlab' filesep 'solvers' filesep 'fuse_sol']);
addpath([rootpath filesep 'matlab' filesep 'simulations']);
addpath([rootpath filesep 'matlab' filesep 'tools']);

% Also include multipol (assume its in the same folder as this repository)
addpath([workspace filesep 'multipol'])

systemsettings.rootpath = rootpath;
