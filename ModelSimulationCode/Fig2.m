% Figure 2. Intertwisted magnitude of noise 
%% define directories
DefineIO;
plot_dir = fullfile(rootdir, 'Prediction','Fig2');
sim_dir = fullfile(rootdir, 'Prediction','Fig2');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end
if ~exist(sim_dir,'dir')
    mkdir(sim_dir);
end
%% loading parrellel CPU cores
Myclust = parcluster();
Npar = Myclust.NumWorkers;
mypool = parpool(Npar);
reps = Npar; % 40; % repetition of simulations to make the results smooth
%% Mixed noise 
V1mean = 88;
V2mean = 83;
V3 = linspace(0, V1mean, 50)';
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
phi = [0, .2, .35, .47, .6, .7, .8, 1]*pi/2;
epsvec = sin(phi)*9;
etavec = cos(phi)*3.63;
products = {'Probability'};
filename = sprintf('Choice_MixedNoise_eps%1.2f_eta%1.2f', max(epsvec), max(etavec));
% simulation
matfile = fullfile(sim_dir, [filename, '.mat']);
if ~exist(matfile, 'file')
    probsa = nan([numel(epsvec), 3, numel(V3)]);
    probsb = nan([numel(epsvec), 3, numel(V3)]);
    for i = 1:numel(epsvec)
        fprintf("Level #%i\n",i);
        eps = epsvec(i);
        eta = etavec(i);
        sdV1 = eps*ones(size(V3));
        sdV2 = eps*ones(size(V3));
        sdV3 = eps*ones(size(V3));
        dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
        pars = [eta, 1, 1, 1];
        tmpa = nan([reps, 3, numel(V3)]);
        tmpb = nan([reps, 3, numel(V3)]);
        parfor ri = 1:reps
            [tmpa(ri,:,:), ~, ~] = dnDNM(dat, pars, 'none', products); % no-constraint model
            [tmpb(ri,:,:), ~, ~] = dnDNM(dat, pars, 'biological', products); % biological model
        end
        probsa(i,:,:) = squeeze(mean(tmpa, 1));
        probsb(i,:,:) = squeeze(mean(tmpb, 1));
    end
    save(matfile, "probsa",  "probsb");
else
    load(matfile);
end

%% visualization
h = figure; hold on;
mycols = [0         0    1.0000
         0    0.3333    0.8333
         0    0.6667    0.6667
         0    1.0000    0.5000
    1.0000    0.7500         0
    1.0000    0.5000         0
    1.0000    0.2500         0
    1.0000         0         0];
ratio = [];
slope = [];
lt = 0.2;
rt = 0.8;
lg = [];
for i = 1:numel(epsvec)
    x = V3/V2mean;
    mask = x >= lt & x <= rt;
    ratio(:,i,1) = squeeze(probsa(i,1,:)./(probsa(i,1,:) + probsa(i,2,:))*100);
    plot(x, ratio(:,i,1), ':', 'Color', mycols(i,:), 'LineWidth', 2);
    ratio(:,i,2) = squeeze(probsb(i,1,:)./(probsb(i,1,:) + probsb(i,2,:))*100);
    lg(i) = plot(x, ratio(:,i,2), '-', 'Color', mycols(i,:), 'LineWidth', 2);
    coefficients = polyfit(x(mask), ratio(mask,i,2), 1);
    slope(i) = coefficients(1);
end
plot([V1mean, V2mean]/V2mean, [1, 1]*min(ratio(:)), 'kv', 'MarkerFaceColor', 'k');
xlim([0, V1mean/V2mean]);
plot([lt, lt], [max(ratio(:)), min(ratio(:))], 'k--');
plot([rt, rt], [max(ratio(:)), min(ratio(:))], 'k--');
mylg = legend(fliplr(lg), {'9.0 0.0',' ',' ',' ',' ',' ',' ','0.0 3.6'}, 'Location', "northeastoutside", 'Box','off');
title(mylg, 'Noise \sigma_{Early} \sigma_{Late}');
xlabel('Scaled V3');
ylabel('% Correct | V1, V2');
mysavefig(h, filename, plot_dir, 12, [5.4, 4]);
%%
h = figure; hold on;
filename = [filename, 'Slope'];
for i = length(slope):-1:1
    bar(9-i, slope(i), 'FaceColor',mycols(i,:));
end
ylabel('Slope');
xlabel('Cases');
mysavefig(h, filename, plot_dir, 12, [2, 2]);

%% Full range panel
filename = sprintf('Choice_MixedNoiseFullrng');
products = {'Probability'};
V1mean = 88;
V2mean = 83;
V3 = linspace(0, V1mean, 50)';
x = V3/V2mean;
mask = x > .2 & x < .8;
epsvec = linspace(0, 13, 101);
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
nsmpls = 1024*1e3;
etavec = linspace(0, 3.63*13/9, 100);
matfile = fullfile(sim_dir, [filename, '.mat']);
if ~exist(matfile, 'file')
    slope = nan([numel(epsvec), numel(etavec)]);
    for j = 1:numel(etavec)
        fprintf("Eta #%i, Loop over Eps: ",j);
        eta = etavec(j);
        for i = 1:numel(epsvec)
            fprintf(".");
            eps = epsvec(i);
            sdV1 = eps*ones(size(V3));
            sdV2 = eps*ones(size(V3));
            sdV3 = eps*ones(size(V3));
            dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
            pars = [eta,eta,eta, 1, 1, 1];
            reps = 40;
            tmp = nan([reps, 3, numel(V3)]);
            parfor ri = 1:reps
                [tmp(ri,:,:), ~, ~] = dnDNM(dat, pars, 'biological', products); % biological model
            end
            probs = squeeze(mean(tmp, 1));
            ratio = probs(:,1)./(probs(:,1) + probs(:,2))*100;
            coefficients = polyfit(x(mask), ratio(mask), 1);
            slope(i,j) = coefficients(1);
        end
        fprintf('\n');
    end
    save(matfile, 'slope');
else
    load(matfile);
end
%%
addpath(genpath('/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/Mylib'));
% please find the redwhiteblue colormap and load it from
% https://www.mathworks.com/matlabcentral/fileexchange/86932-red-white-blue-colormap
h = figure;
subplot(1,2,1); hold on;
imagesc(etavec, epsvec, slope);
phi = [0:0.01:1]*pi/2;
epsline = sin(phi)*9;
etaline = cos(phi)*3.63;
prb = plot(etaline, epsline, '-', 'Color', [85, 85, 85]/255, 'LineWidth', 2);

phis = [0, .2, .35, .47, .6, .7, .8, 1]*pi/2;
epsspots = sin(phis)*9;
etaspots = cos(phis)*3.63;
for ci = 1:numel(phis)
    plot(etaspots(ci), epsspots(ci), '.', 'Color', mycols(ci,:), 'MarkerSize',24);
    plot(etaspots(ci), epsspots(ci), 'o', 'Color', [85, 85, 85]/255, 'MarkerSize', 8);
end

colormap(bluewhitered);
cb = colorbar;
ylabel('\sigma_{Early noise}');
xlabel('\sigma_{Late noise}');
xlim([min(etavec)-.04022, max(etavec)]);
ylim([min(epsvec)-.1, max(epsvec)]);
mysavefig(h, filename, plot_dir, 12, [10, 3.6]);

subplot(1,2,2); hold on;
il = 7;
[etamesh, epsmesh] = meshgrid(etavec, epsvec);
cb = colorbar;

N = 9;
mygray = gray(N);
mag = linspace(.2,2,N);
for magi = 1:N
    epsline = sin(phi)*9*mag(magi);
    etaline = cos(phi)*3.63*mag(magi);
    plot(etaline, epsline, '-', 'Color', mygray(magi,:), 'LineWidth', .5);
end
[Gx, Gy] = gradient(slope);
Gx = RUNAverg(Gx, il, il);
Gy = RUNAverg(Gy, il, il);
norm = sqrt(Gx.^2 + Gy.^2);
rate = (norm.^.4)./(norm);
Gx = Gx.*rate;
Gy = Gy.*rate;
Xd = RUNAverg(etamesh, il, il);
Yd = RUNAverg(epsmesh, il, il);
Ud = Gx;
Vd = Gy;
quiver(Xd,Yd,Ud,Vd,'Color','red');
ylabel('\sigma_{Early noise}');
xlabel('\sigma_{Late noise}');
xlim([min(etavec), max(etavec)]);
ylim([min(epsvec), max(epsvec)]);
mysavefig(h, filename, plot_dir, 12, [9.9, 3.6]);

%%%%%%%%
% The Following code is not included in the current version of the paper,
% Please ignore
%%%%%%%%
%% V3 mean value - variance matrix
V1mean = 88;
V2mean = 83;
eps1 = 4; % early noise for V1
eps2 = 4; % early noise for V2
V3 = linspace(0, V1mean, 100)';
eps3 = linspace(0, 60, 101);
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
sdV1 = eps1*ones(size(V3));
sdV2 = eps2*ones(size(V3));
etavec = [1.0, 1.4286]; % 1, 1.4286; %linspace(.8, 1.9, 8); % different levels of late noise
products = {'Probability'};
for ti = 1:numel(etavec)
    eta = etavec(ti);
    filename = sprintf('Ratio_Model_%iv3max%1.0f_%ivar3max%1.0f_eta%1.2f', numel(V3), max(V3), numel(eps3), max(eps3), eta);
    % simulation
    SimDatafile = fullfile(sim_dir, [filename, '.mat']);
    if ~exist(SimDatafile,'file')
        Ratios = nan(numel(eps3), numel(V3));
        for i = 1:numel(eps3)
            fprintf('Mixed noise - late noise %i/%i, early noise %i/%i\n', ti, numel(etavec), i, numel(eps3));
            sdV3 = eps3(i)*ones(size(V3));
            dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
            pars = [eta, eta, eta, 1, 1, 1];
            tmpb = nan([reps, numel(V3), 3]);
            parfor ri = 1:reps
                [tmpb(ri,:,:), ~, ~] = dnDNM(dat, pars, 'biological', products); % biological model
            end
            probs = squeeze(mean(tmpb, 1));
            Ratios(i,:) = probs(:,1)./(probs(:,1) + probs(:,2))*100;
        end
        save(SimDatafile, "Ratios", "V3", "eta", "V2mean", '-mat');
    end
end
%% visualization
filename = 'CardinalView_Model';
etavec = [1.0, 1.4286];
h = figure;
for ti = 1:2
    eta = etavec(ti);
    simrslts = sprintf('Ratio_Model_%iv3max%1.0f_%ivar3max%1.0f_eta%1.2f', numel(V3), max(V3), numel(eps3), max(eps3), eta);
    SimDatafile = fullfile(sim_dir, [simrslts, '.mat']);
    load(SimDatafile);

    subplot(1,2,ti); hold on;
    colormap("bone");
    cmap = bone(numel(eps3));
    for ri = 1:5:numel(eps3)
        plot(V3/V2mean, Ratios(ri,:), '-', 'Color', cmap(ri,:), 'LineWidth',2);
    end
    plot([V1mean, V2mean]/V2mean, [1, 1]*min(Ratios(:)), 'v', 'MarkerFaceColor', [.7,.7,.7]);
    xlabel('Scaled V3');
    ylabel('% Correct | V1, V2');
    ylim([57, 77]);
    title(sprintf('Late noise eta = %1.2f', eta));
    mysavefig(h, filename, plot_dir, 12, [6, 2.5]);

    % subplot(2,2,2+(ti-1)*2);
    % colormap("jet");
    % imagesc(V3/V2mean, eps3/V2mean, Ratios, [65, 80]);
    % title('');
    % colorbar;
    % xlabel('Scaled V3');
    % ylabel('\sigma_3 (Early noise)');
    % set(gca, 'YDir', 'normal');
    % mysavefig(h, filename, plot_dir, 12, [10, 8]);
end

%% averaging function
function Dp = RUNAverg(D, p, q)
[n, m] = size(D);
ip = 0;
for i = 1:p:(n-p+1)
    ip = ip + 1;
    jp = 0;
    for j = 1:q:(m-q+1)
        jp = jp + 1;
        Dp(ip, jp) = mean(D(i:(i+p-1),j:(j+q-1)),"all");
    end
end
end
