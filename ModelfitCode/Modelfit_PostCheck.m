%% define directories
DefineIO;
%% Loading the data transformed in the code: /Users/bs3667/Noise/modelfit/ModelFit-DataTrnsfrm.m
load(fullfile(Gitdir, 'myData', 'TrnsfrmData.mat'), 'mt');
fitdir = fullfile(rootdir, 'Modelfit');
fit = tdfread(fullfile(fitdir, 'BestRslts.txt'));
plot_dir = fullfile(fitdir, 'plot');
%% Transform data
blacklist = [22102405; 22102705; 22102708; 22071913; 22110306];
sublist = unique(mt.subID);
sublist = sublist(~ismember(sublist, blacklist));
fulllist = unique(mt.subID);
N = length(sublist);
mtconvert = [];
for s = 1:N
    indvtask = mt(mt.subID == sublist(s),:);
    Vtrgt = unique([indvtask.V1; indvtask.V2]);
    mintrgt = min(Vtrgt);
    if mintrgt > 0 % skip this subject if the min target value is zero, because the values cannot be scaled and the value space does not help in testing of the hypothesis
        indvtask.V1scld = indvtask.V1/mintrgt;
        indvtask.V2scld = indvtask.V2/mintrgt;
        indvtask.V3scld = indvtask.V3/mintrgt;
        indvtask.sdV1scld = indvtask.sdV1/mintrgt;
        indvtask.sdV2scld = indvtask.sdV2/mintrgt;
        indvtask.sdV3scld = indvtask.sdV3/mintrgt;
        mtconvert = [mtconvert; indvtask];
    end  
end
mtconvert.choice = mtconvert.chosenItem - 1;
%% Simulation
for modeli = 1:5
    simdat = fullfile(fitdir, sprintf('Model%i_Predict.mat', modeli));
    if ~exist(simdat, 'file')
        mtmodel = [];
        for s = 1:N
            fprintf('Subject %d:\t', s);
            dat = mtconvert(mtconvert.subID == sublist(s), :);
            subjmask = fit.subID == sublist(s);
            Mp = fit.Mp(subjmask & fit.modeli == modeli);
            delta = fit.delta(subjmask & fit.modeli == modeli);
            switch modeli
                case 1
                    %%
                    x = [Mp, delta];
                    probs = McFadden(x, dat);
                    name = 'McFadden';
                case 2
                    scl = fit.scl(subjmask & fit.modeli == modeli);
                    x = [Mp, delta, scl];
                    probs = Mdl2(x, dat);
                    name = 'LinearDistrb';
                case 3
                    wp = fit.wp(subjmask & fit.modeli == modeli);
                    x = [Mp, delta, wp];
                    probs = DN(x, dat);
                    name = 'DN'; %, cut input, independent';
                case 4
                    scl = fit.scl(subjmask & fit.modeli == modeli);
                    wp = fit.wp(subjmask & fit.modeli == modeli);
                    x = [Mp, delta, wp, scl];
                    probs = dDNb(x, dat, 'absorb');
                    name = 'dDNb'; %, cut input, independent';
                case 5
                    scl = fit.scl(subjmask & fit.modeli == modeli);
                    wp = fit.wp(subjmask & fit.modeli == modeli);
                    x = [Mp, delta, wp, scl];
                    probs = dDNd(x, dat, 'absorb');
                    name = 'dDNd'; %, cut SIGMA, independent';
            end
            dat.modelprob1 = gather(probs(:,1));
            dat.modelprob2 = gather(probs(:,2));
            dat.modelprob3 = gather(probs(:,3));
            mtmodel = [mtmodel; dat];
            fprintf('\n');
        end
        mtmodel.ratio = mtmodel.modelprob2./(mtmodel.modelprob1 + mtmodel.modelprob2);
        save(simdat, "mtmodel", '-mat');
        writetable(mtmodel, fullfile(fitdir, sprintf('Model%i_Predict.txt', modeli)), 'Delimiter', '\t');
    else
        load(simdat);
    end
    %% Visualize in sliding windows
    dat = mtmodel(mtmodel.chosenItem ~= 3 & ~isnan(mtmodel.chosenItem),:);
    GrpMean = grpstats(dat, ["subID","TimeConstraint", "Vaguenesscode", "ID3"], "mean", "DataVars", ["V3", "sdV3", "V3scld", "sdV3scld", "choice","ratio"]);
    colorpalette ={'r','#FFBF00','#00FF80','b'};
    rgbMatrix = [
        0, 0, 255;   % Blue
        255, 192, 203; % Pink
        173, 216, 230; % Light Blue
        255, 0, 0     % Red
        ]/255;
    Window = 0.15;
    LowestV3 = 0; %0.2;
    HighestV3 = 1; %.8;
    h = figure;
    filename = sprintf('Model%i_Predict[Full]', modeli);
    vi = 0;
    i = 0;
    for v = [1, 0] % vague, precise
        vi = vi + 1;
        subplot(1,2,vi); hold on;
        for t = [10, 1.5] % low, high
            i = i + 1;
            Ntrial = [];
            choice = [];
            choicese = [];
            ratio = [];
            ratiose = [];
            sdV3scld = [];
            v3vec = LowestV3:.015:HighestV3;
            dat = GrpMean(GrpMean.TimeConstraint == t & GrpMean.Vaguenesscode == v & GrpMean.mean_V3scld >= LowestV3 &  GrpMean.mean_V3scld <= HighestV3,:);
            for v3 = v3vec
                section = dat(dat.mean_V3scld >= v3 - Window & dat.mean_V3scld <= v3 + Window,:);
                Ntrial = [Ntrial, sum(section.GroupCount)];
                choice = [choice, mean(section.mean_choice)];
                choicese = [choicese, std(section.mean_choice)/sqrt(length(section.mean_choice))];
                ratio = [ratio, mean(section.mean_ratio)];
                ratiose = [ratiose, std(section.mean_ratio)/sqrt(length(section.mean_ratio))];
                sdV3scld = [sdV3scld, mean(section.mean_sdV3scld)];
            end
            cut = Ntrial > 400;
            % scatter(v3vec(cut), ratio(cut), Ntrial(cut)/80*5, 'color', colorpalette{i});
            plot(v3vec(cut), choice(cut), '-', 'Color', colorpalette{i}, 'LineWidth', 2);
            plot(v3vec(cut), ratio(cut), 'k--', 'LineWidth', 2);
            % fill([v3vec fliplr(v3vec)], [ratio-ratiose fliplr(ratio+ratiose)], rgbMatrix(vi,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        xlim([LowestV3, HighestV3]);
        xlabel('Scaled V3');
        ylabel('% Correct | V1, V2');
        mysavefig(h, filename, plot_dir, 12, [8, 4]);
    end
    %% Visualization in heatmap
    dat = mtmodel(mtmodel.chosenItem ~= 3 & ~isnan(mtmodel.chosenItem),:);
    GrpMean = grpstats(dat, ["subID", "TimeConstraint", "Vaguenesscode", "ID3"], "mean", "DataVars", ["V3", "sdV3", "V3scld", "sdV3scld", "choice", "ratio"]);
    Window = 0.15;
    Varrng = [min(GrpMean.mean_sdV3scld), .4];% max(GrpMean.mean_sdV3scld)];
    Bindow = 0.15/2;
    h = figure;
    filename = sprintf('ModelFit_%i_Heatmap', modeli);
    ti = 0;
    TimePressure = {'Low','High'};
    for t = [10, 1.5] % low, high
        ti = ti + 1;
        dat = GrpMean(GrpMean.TimeConstraint == t,:);
        v3vec = LowestV3:.03:HighestV3;
        varvec = Varrng(1):.015:Varrng(2);
        Ntrial = NaN(numel(varvec), numel(v3vec));
        ratio = NaN(numel(varvec), numel(v3vec));
        ratiose = NaN(numel(varvec), numel(v3vec));
        sdV3scld = NaN(numel(varvec), numel(v3vec));
        for vi = 1:numel(v3vec)
            for ri = 1:numel(varvec)
                v3 = v3vec(vi);
                r = varvec(ri);
                maskv3 = dat.mean_V3scld >= v3 - Window & dat.mean_V3scld <= v3 + Window;
                maskr3 = dat.mean_sdV3scld >= r - Bindow & dat.mean_sdV3scld <= r + Bindow;
                section = dat(maskv3 & maskr3,:);
                Ntrial(ri,vi) = sum(section.GroupCount);
                ratio(ri,vi) = mean(section.mean_ratio);
                ratiose(ri,vi) = std(section.mean_ratio)/sqrt(length(section.mean_ratio));
                sdV3scld(ri,vi) = mean(section.mean_sdV3scld);
            end
        end
        ratio(Ntrial<50) = NaN;
        subplot(2, 2, 1+(ti-1)*2); hold on;
        colormap("bone");
        cmap = bone(numel(varvec));
        for ri = 1:numel(varvec)
            plot(v3vec, ratio(ri,:), '.-', 'Color', cmap(ri,:));
        end
        title(TimePressure{ti});
        xlabel('Scaled V3');
        ylabel('% Correct | V1 & V2');
        mysavefig(h, filename, plot_dir, 12, [9, 8]);

        subplot(2, 2, 2+(ti-1)*2); hold on;
        colormap("jet");
        imagesc(v3vec, varvec, ratio);
        c = colorbar('Location', 'northoutside');
        % ylim([0,1]);
        ylabel(c, '% Correct | V1 & V2');
        xlabel('Scaled V3');
        ylabel('V3 Variance');
        mysavefig(h, filename, plot_dir, 12, [9, 8]);
    end
end

%% 
function probs = McFadden(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); % late noise difference between time-pressure conditions
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
samples = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
    else
        values = data.(['V',num2str(ci)])';
    end
    samples(:,:,ci) = repmat(values, num_samples, 1);
end
if gpuparallel
    SVs = samples/Mp + gpuArray.randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = gpuArray(data.chosenItem');
else
    SVs = samples/Mp + randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = data.chosenItem';
end
max_from_each_distribution = SVs == max(SVs, [], 3);
probs = squeeze(sum(max_from_each_distribution, 1) / size(SVs, 1));
% nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
% if gpuparallel
%     nll = gather(nll);
% end
end

function probs = Mdl2(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); % late noise difference between time-pressure conditions
scl = x(3); % scaling parameter on early noise
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
Ntrl = size(dat,1);
samples = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
        stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
        samples(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    else
        values = data.(['V',num2str(ci)])';
        stds = data.(['sdV', num2str(ci)])'*scl;
        samples(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    end
end
if gpuparallel
    SVs = samples/Mp + gpuArray.randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = gpuArray(data.chosenItem');
else
    SVs = samples/Mp + randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = data.chosenItem';
end
max_from_each_distribution = SVs == max(SVs, [], 3);
probs = squeeze(sum(max_from_each_distribution, 1) / size(SVs, 1));
% nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
% if gpuparallel
%     nll = gather(nll);
% end
end

function probs = DN(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); % late noise difference between time-pressure conditions
wp = x(3); % do the same transformation on w
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
samples = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
    else
        values = data.(['V',num2str(ci)])';
    end
    samples(:,:,ci) = repmat(values, num_samples, 1);
end
% The product of divisive normalization before adding late noise
DNP = samples./(sum(samples, 3)*wp + Mp);
% adding late noise
if gpuparallel
    SVs = DNP + gpuArray.randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = gpuArray(data.chosenItem');
else
    SVs = DNP + randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = data.chosenItem';
end
max_from_each_distribution = SVs == max(SVs, [], 3);
probs = squeeze(sum(max_from_each_distribution, 1) / size(SVs, 1));
% nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
% if gpuparallel
%     nll = gather(nll);
% end
end

function probs = dDNb(x, dat, mode) % cut inputs, independent
% set the lower boundary for every input value distribution as zero
% samples in the denominator are independent from the numerator
% the SIGMA term in the denominator is non-negative after that.
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); % late noise difference between time-pressure conditions
wp = x(3); % do the same transformation on w
scl = x(4); % scaling parameter on the early noise
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
Ntrl = size(dat,1);
if strcmp(mode, 'absorb')
    samples = [];
    for ci = 1:3
        if gpuparallel
            values = gpuArray(data.(['V',num2str(ci)])');
            stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
            samples(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        else
            values = data.(['V',num2str(ci)])';
            stds = data.(['sdV', num2str(ci)])'*scl;
            samples(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        end
    end
    D1 = [];
    D2 = [];
    D3 = [];
    for ci = 1:3
        if gpuparallel
            values = gpuArray(data.(['V',num2str(ci)])');
            stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
            D1(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D2(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D3(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        else
            values = data.(['V',num2str(ci)])';
            stds = data.(['sdV', num2str(ci)])'*scl;
            D1(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D2(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D3(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        end
    end
elseif strcmp(mode, 'cutoff')
    error('The cutoff boundary algorithm has not been developped yet.');
end
D1 = sum(D1, 3)*wp + Mp;
D2 = sum(D2, 3)*wp + Mp;
D3 = sum(D3, 3)*wp + Mp;
% The product of divisive normalization before adding late noise
DNP = samples./cat(3, D1, D2, D3);
clear D1 D2 D3;
if gpuparallel
    SVs = DNP + gpuArray.randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = gpuArray(data.chosenItem');
else
    SVs = DNP + randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = data.chosenItem';
end
max_from_each_distribution = SVs == max(SVs, [], 3);
probs = squeeze(sum(max_from_each_distribution, 1) / size(SVs, 1));
% nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
% if gpuparallel
%     nll = gather(nll);
% end
end

function probs = dDNd(x, dat, mode) % cut SIGMA, independent
% set the lower boundary for the summed SIGMA in the denominator
% but not for the input values in the numerator
% samples in the denominator are independent from the numerator
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); % late noise difference between time-pressure conditions
wp = x(3); % do the same transformation on w
scl = x(4); % scaling parameter on the early noise
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
Ntrl = size(dat,1);
samples = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
        stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
        samples(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    else
        values = data.(['V',num2str(ci)])';
        stds = data.(['sdV', num2str(ci)])'*scl;
        samples(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    end
end
D1 = [];
D2 = [];
D3 = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
        stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
        D1(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        D2(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        D3(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    else
        values = data.(['V',num2str(ci)])';
        stds = data.(['sdV', num2str(ci)])'*scl;
        D1(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        D2(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        D3(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
    end
end
if strcmp(mode, 'absorb')
    D1 = max(sum(D1, 3),0)*wp + Mp;
    D2 = max(sum(D2, 3),0)*wp + Mp;
    D3 = max(sum(D3, 3),0)*wp + Mp;
elseif strcmp(mode, 'cutoff')
    error('The cutoff boundary algorithm has not been developped yet.');
end
% The product of divisive normalization before adding late noise
DNP = samples./cat(3, D1, D2, D3);
clear D1 D2 D3;
if gpuparallel
    SVs = DNP + gpuArray.randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = gpuArray(data.chosenItem');
else
    SVs = DNP + randn(size(samples)).*(1 + delta*repmat(data.TimeConstraint'==1.5,num_samples,1,3))*eta;
    choice = data.chosenItem';
end
max_from_each_distribution = SVs == max(SVs, [], 3);
probs = squeeze(sum(max_from_each_distribution, 1) / size(SVs, 1));
% nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
% if gpuparallel
%     nll = gather(nll);
% end
end
