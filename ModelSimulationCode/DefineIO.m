%% define directories
[os, ~, ~] = computer;
if strcmp(os,'MACA64')
    rootdir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/NoiseProject';
    Gitdir = '~/Noise';
elseif strcmp(os,'GLNXA64')
    rootdir = '/gpfs/data/glimcherlab/BoShen/NoiseProject';
    % rootdir = '/scratch/bs3667/NoiseProject';
    Gitdir = '~/Noise';
elseif strcmp(os, 'PCWIN64')
    rootdir = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\NoiseProject';
    Gitdir = 'C:\Users\Bo\Documents\GitHub\Noise';
end
addpath(genpath(Gitdir));
