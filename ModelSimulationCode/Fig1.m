% Figure 1. Generating the predictions for the early and late noise
%% define directories
DefineIO;
plot_dir = fullfile(rootdir, 'Prediction','Fig1');
sim_dir = fullfile(rootdir, 'Prediction','Fig1');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end
if ~exist(sim_dir,'dir')
    mkdir(sim_dir);
end

%% loading parrellel CPU cores
Myclust = parcluster();
Npar = Myclust.NumWorkers;
mypool = parpool(20);
reps = 40; % repetition of simulations to make the results smooth
%% Early and late noise only
filename = sprintf('SNR_Ovlp_Choice_v3');
V1mean = 88;
V2mean = 83;
V3 = linspace(0 , V1mean, 50)';
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
products = {'Probability','Coeff_of_Var','Overlap'};
% simulation
matfile = fullfile(sim_dir, [filename, '.mat']);
if ~exist(matfile, 'file')
    %% with early noise only
    eps = 9;
    eta = 0;
    sdV1 = eps*ones(size(V3));
    sdV2 = eps*ones(size(V3));
    sdV3 = eps*ones(size(V3));
    dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
    pars = [eta,eta,eta, 1, 1, 1];
    tmp1a = nan([reps, numel(V3), 3]); % probs
    tmp2a = nan([reps, numel(V3)]); % Ovlps
    tmp3a = nan([reps, numel(V3), 3]); % CVs
    tmp1b = nan([reps, numel(V3), 3]); % probs
    tmp2b = nan([reps, numel(V3)]); % Ovlps
    tmp3b = nan([reps, numel(V3), 3]); % CVs
    parfor ri = 1:reps
        fprintf('Early noise, rep %i\n', ri);
        [tmp1a(ri,:,:), tmp2a(ri,:), tmp3a(ri,:,:)] = dnDNM(dat, pars, 'none', products); % no-constraint model
        fprintf('.');
        [tmp1b(ri,:,:), tmp2b(ri,:), tmp3b(ri,:,:)] = dnDNM(dat, pars, 'biological', products); % biological model
        fprintf('.\n');
    end
    probsa = squeeze(mean(tmp1a, 1));
    Ovlpsa = squeeze(mean(tmp2a, 1));
    CVsa = squeeze(mean(tmp3a, 1));
    probsb = squeeze(mean(tmp1b, 1));
    Ovlpsb = squeeze(mean(tmp2b, 1));
    CVsb = squeeze(mean(tmp3b, 1));

    %% with late noise only
    eps = 0;
    eta = 4;
    sdV1 = eps*ones(size(V3));
    sdV2 = eps*ones(size(V3));
    sdV3 = eps*ones(size(V3));
    dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
    pars = [eta,eta,eta, 1, 1, 1];
    tmp1a = nan([reps, numel(V3),3]);
    tmp2a = nan([reps, numel(V3)]);
    tmp3a = nan([reps, numel(V3), 3]);
    tmp1b = nan([reps, numel(V3), 3]);
    tmp2b= nan([reps, numel(V3)]);
    tmp3b = nan([reps, numel(V3), 3]);
    parfor ri = 1:reps
        fprintf('Late noise, rep %i\n', ri);
        [tmp1a(ri,:,:), tmp2a(ri,:), tmp3a(ri,:,:)] = dnDNM(dat, pars, 'none', products); % no-constraint model
        fprintf('.');
        [tmp1b(ri,:,:), tmp2b(ri,:), tmp3b(ri,:,:)] = dnDNM(dat, pars, 'biological', products); % biological model
        fprintf('.\n');
        % output: probs 50*3, Ovlps 1*50, CVs 50*3
    end
    probsaL = squeeze(mean(tmp1a, 1));
    OvlpsaL = squeeze(mean(tmp2a, 1));
    CVsaL = squeeze(mean(tmp3a, 1));
    probsbL = squeeze(mean(tmp1b, 1));
    OvlpsbL = squeeze(mean(tmp2b, 1));
    CVsbL = squeeze(mean(tmp3b, 1));
    save(matfile, "probsa",  "probsb", "Ovlpsa", "Ovlpsb", "CVsa", "CVsb",  "probsaL", "probsbL", "OvlpsaL", "OvlpsbL", "CVsaL", "CVsbL");
else
    load(matfile);
end
%% plotting
h = figure;
subplot(3,2,1); hold on;
plot(V3, 1./CVsa(:,1), 'r:', 'LineWidth', 2);
plot(V3, 1./CVsb(:,1), 'r-', 'LineWidth', 2);
plot(V3, 1./CVsa(:,2), 'r:', 'LineWidth', 2);
plot(V3, 1./CVsb(:,2), 'r-', 'LineWidth', 2);
xlabel('V3');
ylabel('Single-item SNR');
xlim([0, V1mean]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%0.3f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
legend({'Theoritical','Biological'}, 'Location', 'northeast', 'FontSize',10);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

subplot(3,2,2); hold on;
plot(V3, 1./CVsaL(:,1), 'c:', 'LineWidth', 2);
plot(V3, 1./CVsbL(:,1), 'c-', 'LineWidth', 2);
plot(V3, 1./CVsaL(:,2), 'c:', 'LineWidth', 2);
plot(V3, 1./CVsbL(:,2), 'c-', 'LineWidth', 2);
xlabel('V3');
ylabel('Single-item SNR');
xlim([0, V1mean]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%0.3f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

subplot(3,2,3); hold on;
plot(V3, Ovlpsa, 'r:', 'LineWidth', 2);
plot(V3, Ovlpsb, 'r-', 'LineWidth', 2);
xlabel('V3');
ylabel('% Overlap | V1, V2');
xlim([0, V1mean]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%2.1f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

subplot(3,2,4); hold on;
plot(V3, OvlpsaL, 'c:', 'LineWidth', 2);
plot(V3, OvlpsbL, 'c-', 'LineWidth', 2);
xlabel('V3');
ylabel('% Overlap | V1, V2');
xlim([0, V1mean]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%2.1f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

subplot(3,2,5); hold on;
ratio = probsa(:,1)./(probsa(:,1) + probsa(:,2))*100;
plot(V3, ratio, 'r:', 'LineWidth', 2);
plot([V1mean, V2mean], [1, 1]*min(ratio), 'kv', 'MarkerFaceColor', [.7, .7, .7]);
ratio = probsb(:,1)./(probsb(:,1) + probsb(:,2))*100;
plot(V3, ratio, 'r-', 'LineWidth', 2);
xlabel('V3');
ylabel('% Correct | V1, V2');
xlim([0, V1mean]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%2.1f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

subplot(3,2,6); hold on;
ratio = probsaL(:,1)./(probsaL(:,1) + probsaL(:,2))*100;
plot(V3, ratio, 'c:', 'LineWidth', 2);
plot([V1mean, V2mean], [1, 1]*min(ratio), 'kv', 'MarkerFaceColor', [.7, .7, .7]);
ylim([61, 65]);
ratio = probsbL(:,1)./(probsbL(:,1) + probsbL(:,2))*100;
plot(V3, ratio, 'c-', 'LineWidth', 2);
xlabel('V3');
ylabel('% Correct | V1, V2');
xlim([0, V1mean]);
% ylim([61,65]);
% yticks = get(gca, 'YTick');
% yticklabels = arrayfun(@(v) sprintf('%2.1f', v), yticks, 'UniformOutput', false);
% set(gca, 'YTickLabel', yticklabels);
mysavefig(h, filename, plot_dir, 12, [6, 8]);

%% Test a cardinal view of V3 magnitude - V3 variance 
% filename = sprintf('V3mag100_var101_Choice_Ovlp');
filename = sprintf('V3mag50_var6_Choice_Ovlp_II');
matfile = fullfile(sim_dir, [filename, '.mat']);
products = {'Probability','Overlap'}; % 'Coeff_of_Var',
V1mean = 88;
V2mean = 83;
% V3 = linspace(0, V1mean, 100)';
V3 = linspace(0 , V1mean, 50)';
% eps3 = linspace(0, 18, 101);
eps3 = linspace(0, 14, 6);
V1 = V1mean*ones(size(V3));
eps1 = 9;
V2 = V2mean*ones(size(V3));
eps2 = 9;
sdV1 = eps1*ones(size(V3));
sdV2 = eps2*ones(size(V3));
eta = 0;
% simulation
if ~exist(matfile, 'file')
    probsa = nan([numel(eps3), numel(V3), 3]);
    Ovlpsa = nan([numel(eps3), numel(V3)]);
    probsb = nan([numel(eps3), numel(V3), 3]);
    Ovlpsb = nan([numel(eps3), numel(V3)]);
    for i = 1:numel(eps3)
        fprintf('eps3 %i/%i\n', i, numel(eps3));
        sdV3 = eps3(i)*ones(size(V3));
        dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
        pars = [eta,eta,eta, 1, 1, 1];
        tmp1a = nan([reps, numel(V3), 3]);
        tmp2a = nan([reps, numel(V3)]);
        tmp1b = nan([reps, numel(V3), 3]);
        tmp2b = nan([reps, numel(V3)]);
        parfor ri = 1:reps
            fprintf('Early noise, rep %i\n', ri);
            [tmp1a(ri,:,:), tmp2a(ri,:), ~] = dnDNM(dat, pars, 'none', products); % no-constraint model
            [tmp1b(ri,:,:), tmp2b(ri,:), ~] = dnDNM(dat, pars, 'biological', products); % biological model
            % output: probs, Ovlps, CVs
        end
        probsa(i,:,:) = squeeze(mean(tmp1a, 1));
        Ovlpsa(i,:) = squeeze(mean(tmp2a, 1));
        probsb(i,:,:) = squeeze(mean(tmp1b, 1));
        Ovlpsb(i,:) = squeeze(mean(tmp2b, 1));
    end
    save(matfile, "probsa", "Ovlpsa", "probsb", "Ovlpsb");
else
    load(matfile);
end
%% plotting 
h = figure;
figsz = [3, 4.8];
exmpls = 1:6;
ncols = round(numel(exmpls)*1.2);
reds = flip([ones(ncols,1), linspace(0,1,ncols)', linspace(0,1,ncols)']); % Gradient from red to white
subplot(2,1,1); hold on;
for ri = 1:numel(exmpls)
    plot(V3, Ovlpsa(exmpls(ri),:), ':', "Color", reds(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
    plot(V3, Ovlpsb(exmpls(ri),:), '-', "Color", reds(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
end
plot([V1mean, V2mean], [1, 1]*min(Ovlpsb(:)), 'v', 'MarkerFaceColor', [.7,.7,.7]);
ylabel('% Overlap | V1, V2');
xlabel('V3');
xlim([min(V3), max(V3)]);
mysavefig(h, filename, plot_dir, 12, figsz);

ratioa = probsa(:,:,1)./(probsa(:,:,1) + probsa(:,:,2))*100;
ratiob = probsb(:,:,1)./(probsb(:,:,1) + probsb(:,:,2))*100;
x = V3/V2mean;
lt = 0.2;
rt = 0.8;
mask = x >= lt & x <= rt;
slope = [];
lgd = [];
subplot(2,1,2); hold on;
for ri = 1:numel(exmpls)
    plot(V3, ratioa(exmpls(ri),:), ':', "Color", reds(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
    lgd(ri) = plot(V3, ratiob(exmpls(ri),:), '-', "Color", reds(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
    coefficients = polyfit(x(mask), ratiob(exmpls(ri),mask), 1);
    slope(ri) = coefficients(1);
end
plot(lt*V2mean*[1, 1], [64, 65], 'k--', 'LineWidth',.5);
plot(rt*V2mean*[1, 1], [64, 65], 'k--', 'LineWidth',.5);
plot([V1mean, V2mean], [1, 1]*61, 'v', 'MarkerFaceColor', [.7,.7,.7]);
numericVector = eps3(exmpls);
cellArray = arrayfun(@(x) sprintf('%.1f', x), numericVector, 'UniformOutput', false);
legend(lgd, cellArray, 'Location', 'best');
ylabel('% Correct | V1, V2');
xlabel('V3');
xlim([min(V3), max(V3)]);
mysavefig(h, filename, plot_dir, 12, figsz);
%%
h = figure; hold on;
filename = [filename, 'Slope'];
for ri = 1:numel(exmpls)
    bar(ri, slope(ri), 'FaceColor', reds(ri + (ncols-numel(exmpls)),:));
end
ylabel('Slope');
xlabel('Cases');
mysavefig(h, filename, plot_dir, 12, [2, 2]);
%% Test a cardinal view of V3 magnitude - V3 lalte noise  
filename = sprintf('V3mag50_latenoise6_Choice_Ovlp');
matfile = fullfile(sim_dir, [filename, '.mat']);
products = {'Probability','Overlap'}; % 'Coeff_of_Var',
V1mean = 88;
V2mean = 83;
V3 = linspace(0, V1mean, 50)';
eta3s = linspace(0, 6, 6)';
V1 = V1mean*ones(size(V3));
eta1 = 4;
V2 = V2mean*ones(size(V3));
eta2 = 4;
eps1 = 0;
eps2 = 0;
eps3 = 0;
sdV1 = eps1*ones(size(V3));
sdV2 = eps2*ones(size(V3));
sdV3 = eps3*ones(size(V3));
% simulation
if ~exist(matfile, 'file')
    probsb = nan([numel(eta3s), numel(V3), 3]);
    Ovlpsb = nan([numel(eta3s), numel(V3)]);
    for i = 1:numel(eta3s)
        fprintf('eta3 %i/%i\n', i, numel(eta3s));
        dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
        pars = [eta1,eta2,eta3s(i), 1, 1, 1];
        tmp2b = nan([reps, numel(V3)]);
        parfor ri = 1:reps
            fprintf('Late noise, rep %i\n', ri);
            [tmp1b(ri,:,:), tmp2b(ri,:), ~] = dnDNM(dat, pars, 'biological', products); % biological model
            % output: probs, Ovlps, CVs
        end
        probsb(i,:,:) = squeeze(mean(tmp1b, 1));
        Ovlpsb(i,:) = squeeze(mean(tmp2b, 1));
    end
    save(matfile, "probsb", "Ovlpsb");
else
    load(matfile);
end
%% plotting 
h = figure;
exmpls = 1:6;
ncols = round(numel(exmpls)*1.2);
blues = flip([linspace(0,1,ncols)', linspace(0,1,ncols)', ones(ncols,1)]); % Gradient from red to white
subplot(2,1,1); hold on;
for ri = 1:numel(exmpls)
    plot(V3, Ovlpsb(exmpls(ri),:), '-', "Color", blues(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
end
plot([V1mean, V2mean], [1, 1]*min(Ovlpsb(:)), 'v', 'MarkerFaceColor', [.7,.7,.7]);
ylabel('% Overlap | V1, V2');
xlabel('V3');
xlim([min(V3), max(V3)]);
mysavefig(h, filename, plot_dir, 12, figsz);

ratiob = probsb(:,:,1)./(probsb(:,:,1) + probsb(:,:,2))*100;
x = V3/V2mean;
lt = 0.2;
rt = 0.8;
mask = x >= lt & x <= rt;
slope = [];
lgd = [];
subplot(2,1,2); hold on;
for ri = 1:numel(exmpls)
    lgd(ri) = plot(V3, ratiob(exmpls(ri),:), '-', "Color", blues(ri + (ncols-numel(exmpls)),:), 'LineWidth', 2);
    coefficients = polyfit(x(mask), ratiob(exmpls(ri),mask), 1);
    slope(ri) = coefficients(1);
end
plot([V1mean, V2mean], [1, 1]*min(ratiob(:)), 'v', 'MarkerFaceColor', [.7,.7,.7]);
numericVector = eta3s(exmpls);
cellArray = arrayfun(@(x) sprintf('%.1f', x), numericVector, 'UniformOutput', false);
legend(lgd, cellArray, 'Location', 'best');
ylabel('% Correct | V1, V2');
xlabel('V3');
xlim([min(V3), max(V3)]);
mysavefig(h, filename, plot_dir, 12, figsz);
%%
h = figure; hold on;
filename = [filename, 'Slope'];
for ri = 1:numel(exmpls)
    bar(ri, slope(ri), 'FaceColor', blues(ri + (ncols-numel(exmpls)),:));
end
ylabel('Slope');
xlabel('Cases');
mysavefig(h, filename, plot_dir, 12, [2, 2]);
