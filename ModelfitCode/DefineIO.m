%% define directories
[os, ~, ~] = computer;
if strcmp(os,'MACA64') % if on local Mac OS
    rootdir = fullfile('..', 'OutputfromMatlab');
    Gitdir = '../';
elseif strcmp(os, 'PCWIN64') % if on local Windows
    rootdir = fullfile('..', 'OutputfromMatlab');
    Gitdir = '..';
elseif strcmp(os,'GLNXA64') % if on remote HPC
    rootdir = '/gpfs/data/glimcherlab/BoShen/NoiseProject';
    Gitdir = '~/Shen2025NatComm';
end

addpath(genpath(Gitdir));
datadir = fullfile(Gitdir,'myData');

Fitdir = fullfile(rootdir, 'Modelfit');
if ~exist(Fitdir, 'dir')
    mkdir(Fitdir);
end
plotdir = fullfile(Fitdir, 'plot');
if ~exist(plotdir, 'dir')
    mkdir(plotdir);
end
mtrxdir = fullfile(Fitdir, 'Objs');
if ~exist(mtrxdir, 'dir')
    mkdir(mtrxdir);
end

