% Figure 4 Extend. Visualize the predictions from alternative models for revision
%% define directories
DefineIO;
plot_dir = fullfile(rootdir, 'Prediction','Fig4');
sim_dir = fullfile(rootdir, 'Prediction','Fig4');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end
if ~exist(sim_dir,'dir')
    mkdir(sim_dir);
end
mycols = {'#00FF80','#FF0000';
    '#0000FF','#FFBF00'}; 
%% load data
load(fullfile(Gitdir,'myData','TrnsfrmData.mat'));
sublist = unique(mt.subID);
blacklist = [22102405; 22102705; 22102708; 22071913; 22110306];
sublist = sublist(~ismember(sublist, blacklist));
N = length(sublist);
mtconvert = [];
for s = 1:N
    indvtask = mt(mt.subID == sublist(s),:);
    Vtrgt = unique([indvtask.V1; indvtask.V2]);
    mintrgt = min(Vtrgt);
    if mintrgt > 0 % skip this subject if the min target value is zero, because the values cannot be scaled and the value space does not help in testing of the hypothesis
        indvtask.V1 = indvtask.V1/mintrgt;
        indvtask.V2 = indvtask.V2/mintrgt;
        indvtask.V3 = indvtask.V3/mintrgt;
        indvtask.sdV1 = indvtask.sdV1/mintrgt;
        indvtask.sdV2 = indvtask.sdV2/mintrgt;
        indvtask.sdV3 = indvtask.sdV3/mintrgt;
        mtconvert = [mtconvert; indvtask];
    end  
end
mtconvert.choice = mtconvert.chosenItem - 1;
%% V3 variance over sliding window of scaled V3
V3scld = linspace(0, 1, 50)';
step = .05;
wndw = .3;
meansdV3 = [];
Vagueness = {'Precise','Vague'};
TimePressure = {'Low','High'};
h = figure; hold on;
for i = 1:2 
    for ti = 1:2
        for v3i = 1:numel(V3scld)
            mask = strcmp(mtconvert.TimePressure,TimePressure{ti}) & ...
                strcmp(mtconvert.Vagueness, Vagueness{i}) & ...
                mtconvert.V3 >= V3scld(v3i)-wndw/2 & ...
                mtconvert.V3 <= V3scld(v3i)+wndw/2;
            meansdV3(v3i,ti,i) = mean(mtconvert.sdV3(mask));
        end
        plot(V3scld, meansdV3(:,ti,i),'Color',mycols{ti,i});
    end
end
%% loading parrellel CPU cores
Myclust = parcluster();
Npar = Myclust.NumWorkers;
mypool = parpool(Npar/2);
reps = 40; % repetition of simulations to make the results smooth

%% graded color, two panels
V1mean = 88;
V2mean = 83;
eps1 = 4.5; % early noise for V1
eps2 = 4.5; % early noise for V2
V3 = linspace(0, V2mean, 50)';
V3mean = mean(V3);
eps3vec = linspace(0, .35, 6)*V2mean;
V1 = V1mean*ones(size(V3));
V2 = V2mean*ones(size(V3));
sdV1 = eps1*ones(size(V3));
sdV2 = eps2*ones(size(V3));
etavec = [1, 1.4286]; % multiple levels of late noise
K = 75;
products = {'Probability'};
for modeli = 1:4
    filename = sprintf('Ratio_Model%i_%iv3max%1.0f_%s', modeli, numel(V3), max(V3), '6lines');
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
            sdV3 = eps3vec(i)*ones(size(V3));
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
