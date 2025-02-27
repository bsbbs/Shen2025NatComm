% Mean-scaled noise. Generating the predictions for the early and late
% noise with mean-scaled noise
%% define directories
DefineIO;
plot_dir = fullfile(rootdir, 'Prediction','MeanScaledNoise');
sim_dir = fullfile(rootdir, 'Prediction','MeanScaledNoise');
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
    sdV3 = eps*V3/V2mean;
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
    save(matfile, "probsa",  "Ovlpsa",  "CVsa", "probsb",  "Ovlpsb",  "CVsb");
else
    load(matfile);
end

%% graded color, two panels
V1mean = 88;
V2mean = 83;
eps1 = 4.5; % early noise for V1
eps2 = 4.5; % early noise for V2
V3 = linspace(0, V2mean, 50)';
V3mean = mean(V3);
eps3vec = linspace(0, .1084, 6);%*V2mean;
eps3vec = sort([eps3vec, 4.5/V2mean]);
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
sdV1 = eps1*ones(size(V3));
sdV2 = eps2*ones(size(V3));
etavec = [0, 1.4286]; % multiple levels of late noise
K = 75;
products = {'Probability'};
for modeli = 1:4
    filename = sprintf('Ratio_Model%i_MeanScaled_V1_%i_sd1_%1.1f_%iV3max%1.0f_%s', modeli, V1mean, eps1, numel(V3), max(V3));
    Rslts = table('Size', [0 4], 'VariableTypes', {'double', 'double', 'double', 'double'},...
    'VariableNames', {'Early', 'Late', 'V3', 'choice'});
    SimDatafile = fullfile(sim_dir, [filename, '.mat']);
    % simulation
    if exist(SimDatafile,'file')
        load(SimDatafile);
    else
        fprintf('Simulate Model %i of 4\n', modeli);
        switch modeli
            case 1 % independent linear coding model, late noise
                % controller = [(eps1+eps2)/2/(1+V1mean+V2mean+V3mean)*K, 0, 0, 1+V1mean+V2mean+V3mean]; 
                controller = [0, 0, 0, 1+V1mean+V2mean+V3mean]; 
                %  1. equivalence of early noise to late noise; 2. scaling for early noise; 3. weight of normalization; 4. baseline normalization
            case 2 % independent linear coding model, early + late noise
                controller = [0, 1, 0, 1+V1mean+V2mean+V3mean];
            case 3 % divisive normalization model, late noise
                % controller = [(eps1+eps2)/2/(1+V1mean+V2mean+V3mean)*K, 0, 1, 1];
                controller = [0, 0, 1, 1];
            case 4 % divisive normalization model, early + late noise
                controller = [0, 1, 1, 1];
        end
        Ratios = nan(numel(eps3vec), numel(etavec), numel(V3));
        for i = 1:numel(eps3vec)
            fprintf('Early noise %i of %i\n', i, numel(eps3vec));
            sdV3 = eps3vec(i)*V3;
            dat = table(V1,V2,V3,sdV1,sdV2,sdV3);
            for ti = 1:numel(etavec)
                fprintf('\tLate noise %i of 2', ti);
                eta = etavec(ti);
                pars = [ones(1,3)*(eta+controller(1)), controller(2:4)];
                tmpb = nan([reps, numel(V3), 3]);
                parfor ri = 1:reps
                    [tmpb(ri,:,:), ~, ~] = dnDNM(dat, pars, 'biological', products); % biological model
                    fprintf('.');
                end
                fprintf('\n');
                probs = squeeze(mean(tmpb, 1));
                Ratios(i,ti,:) = probs(:,1)./(probs(:,1) + probs(:,2))*100;
            end
        end
        xval = V3'/V2mean;
        save(SimDatafile, "Ratios","xval",'-mat');
    end
    for i = 1:numel(eps3vec)
        for ti = 1:numel(etavec)
            new_row = table(repmat(eps3vec(i),numel(V3),1), repmat(etavec(ti),numel(V3),1), V3, squeeze(Ratios(i,ti,:)),'VariableNames', Rslts.Properties.VariableNames);
            Rslts = [Rslts; new_row];
        end
    end
    writetable(Rslts, fullfile(sim_dir, [filename, '.txt']), 'Delimiter', '\t');

    %% visualization
    xval = V3'/V2mean;
    lt = 0.2;
    rt = 0.8;
    mask = xval >= lt & xval <= rt;
    cmap = jet(numel(eps3vec));
    h = figure;
    for ti = 1:numel(etavec)
        if modeli == 4
            subplot(1,2,ti);
        end
        hold on;
        lg = [];
        for i = 1:numel(eps3vec)
            ratio = squeeze(Ratios(i, ti, :))';
            lg(ti,i) = plot(xval, ratio, '-', 'LineWidth', 2, 'Color', cmap(i,:));
        end
        if ti == 1 || modeli == 4
            plot([V1mean, V2mean]/V2mean, [1, 1]*min(ratio(:)), 'kv', 'MarkerFaceColor', [.7,.7,.7], 'MarkerSize', 5);
            plot([lt, lt], [ylim], 'k--');
            plot([rt, rt], [ylim], 'k--');
            xlabel('Scaled V3');
            ylabel('% Correct | V1, V2');
        end
        if modeli == 4
            mysavefig(h, filename, plot_dir, 12, [7.27, 2.37]);
        end
    end
    if modeli < 4
        mysavefig(h, filename, plot_dir, 12, [7.27/2, 2.37]);
    end
end