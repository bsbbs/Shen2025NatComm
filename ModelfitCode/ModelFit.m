%% define directoriress
DefineIO;

% Switch to the working directory
cd(Fitdir);

%% load data
load(fullfile(datadir, 'TrnsfrmData.mat'));
% blacklist = [22102401, 22102405, 22110306]; % these subjects report they aimed to choose smaller-value items.
% mt = mt(~ismember(mt.subID, blacklist),:);
sublist = unique(mt.subID);
blacklist = [22102405; 22102705; 22102708; 22071913; 22110306];
% disp(head(mt));
%% Maximum likelihood fitting to the choice behavior
options = bads('defaults');     % Default options
options.Display = 'final';
options.UncertaintyHandling = true;    %s Function is stochastic
options.NoiseFinalSamples = 30;
Rslts = table('Size', [0 12], 'VariableTypes', {'double', 'double', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'logical', 'uint16'},...
    'VariableNames', {'subID', 'modeli', 'name', 'Mp', 'delta', 'wp', 'scl', 'delta2', 'nll', 'nllsd', 'success', 'iterations'});
testfile = fullfile(Fitdir, 'AllRslts.txt');
if ~exist(testfile, 'file')
    fp = fopen(testfile, 'w+');
    fprintf(fp, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subID', 'Model', 'randi', 'Mp0', 'delta_0', 'wp0', 'scl0', 'detla2_0', 'Mp', 'delta', 'wp', 'scl', 'detla2', 'nll', 'nllsd', 'success', 'iterations');
    fclose(fp);
end
Myclust = parcluster();
Npar = Myclust.NumWorkers;
fitlist = find(~ismember(sublist, blacklist));
for subj = fitlist'
    %%
    fprintf('Subject %d:\n', subj);
    dat = mt(mt.subID == sublist(subj), :);
    for modeli = 1:5
        switch modeli
            case 1
                nLLfunc = @(x) McFadden(x, dat);
                name = 'McFadden';
            case 2
                nLLfunc = @(x) Mdl2(x, dat);
                name = 'dnLinear';
            case 3
                nLLfunc = @(x) DNM(x, dat);
                name = 'DNM'; %, cut input, independent';
            case 4
                nLLfunc = @(x) dnDNM(x, dat, 'absorb');
                name = 'dnDNM'; %, cut input, independent';
            case 5
                nLLfunc = @(x) dnDNMe(x, dat, 'absorb');
                name = 'dnDNMExtened'; %, cut input, independent';
        end
        fprintf('\tModel %i\n', modeli);
        filename = fullfile(mtrxdir, sprintf('Subj%02i_Mdl%i.mat', subj, modeli));
        if ~exist(filename, 'file')
            if isempty(gcp("nocreate"))
                mypool = parpool(Npar);
            end
            % shared parameters for all models
            LB = [0, -1]; % [Mp, delta1]
            UB = [1000, 4];
            PLB = [1, -.2];
            PUB = [100, 2];
            % nest other parameters
            if modeli == 2 % nest scaling on early noise for model 2
                LB = [LB, 0]; % [Mp, delta1, scl]
                UB = [UB, 4];
                PLB = [PLB, 0];
                PUB = [PUB, 2];
            end
            if modeli == 3 % nest divisive normalization for model 3
                LB = [LB, 0]; % [Mp, delta1, wp]
                UB = [UB, 4];
                PLB = [PLB, 0];
                PUB = [PUB, 1.4];
            end
            if modeli == 4 % nest for model 4
                LB = [LB, 0, 0]; % [Mp, log(1+delta), wp, scl]
                UB = [UB, 4, 4];
                PLB = [PLB, 0, 0];
                PUB = [PUB, 1.4, 2];
            end
            if modeli == 5 % nest for model 5
                LB = [LB, 0, 0, -1]; % [Mp, delta1, wp, scl, delta2]
                UB = [UB, 4, 4, 4];
                PLB = [PLB, 0, 0, -.2];
                PUB = [PUB, 1.4, 2, 2];
            end

            nLL = [];
            params = {};
            success = [];
            res = {};
            parfor i = 1:Npar
                x0 = PLB + (PUB - PLB) .* rand(size(PLB));
                [xOpt,fval,exitflag,output] = bads(nLLfunc,x0,LB,UB,PLB,PUB,[],options);
                % 'Mp', 'delta', 'wp', 'scl',
                if modeli == 1
                    dlmwrite(testfile, [sublist(subj), modeli, i, x0, NaN, NaN, NaN, xOpt, NaN, NaN, NaN, fval, output.fsd, exitflag, output.iterations],'delimiter','\t','precision','%.6f','-append');
                elseif modeli == 2
                    dlmwrite(testfile, [sublist(subj), modeli, i, x0(1:2), NaN, x0(3), NaN, xOpt(1:2), NaN, xOpt(3), NaN, fval, output.fsd, exitflag, output.iterations],'delimiter','\t','precision','%.6f','-append');
                elseif modeli == 3
                    dlmwrite(testfile, [sublist(subj), modeli, i, x0, NaN, NaN, xOpt, NaN, NaN, fval, output.fsd, exitflag, output.iterations],'delimiter','\t','precision','%.6f','-append');
                elseif modeli == 4
                    dlmwrite(testfile, [sublist(subj), modeli, i, x0, NaN, xOpt, NaN, fval, output.fsd, exitflag, output.iterations],'delimiter','\t','precision','%.6f','-append');
                elseif modeli == 5
                    dlmwrite(testfile, [sublist(subj), modeli, i, x0, xOpt, fval, output.fsd, exitflag, output.iterations],'delimiter','\t','precision','%.6f','-append');
                end
                params{i} = xOpt;
                nLL(i) = fval;
                success(i) = exitflag;
                res{i} = output;
            end
            besti = find(nLL == min(nLL));
            xOpt = params{besti};
            fval = nLL(besti);
            exitflag = success(besti);
            output = res{besti};
            save(filename, 'xOpt', 'fval', 'exitflag', 'output');
        else
            load(filename);
        end
        if modeli == 1
            new_row = table(sublist(subj), modeli, {name}, xOpt(1), xOpt(2), NaN, NaN, NaN, fval, output.fsd, exitflag, output.iterations, 'VariableNames', Rslts.Properties.VariableNames);
        elseif modeli == 2
            new_row = table(sublist(subj), modeli, {name}, xOpt(1), xOpt(2), NaN, xOpt(3), NaN, fval, output.fsd, exitflag, output.iterations, 'VariableNames', Rslts.Properties.VariableNames);
        elseif modeli == 3
            new_row = table(sublist(subj), modeli, {name}, xOpt(1), xOpt(2), xOpt(3), NaN, NaN, fval, output.fsd, exitflag, output.iterations, 'VariableNames', Rslts.Properties.VariableNames);
        elseif modeli == 4
            new_row = table(sublist(subj), modeli, {name}, xOpt(1), xOpt(2), xOpt(3), xOpt(4), NaN, fval, output.fsd, exitflag, output.iterations, 'VariableNames', Rslts.Properties.VariableNames);
        elseif modeli == 5
            new_row = table(sublist(subj), modeli, {name}, xOpt(1), xOpt(2), xOpt(3), xOpt(4), xOpt(5), fval, output.fsd, exitflag, output.iterations, 'VariableNames', Rslts.Properties.VariableNames);
        end
        Rslts = [Rslts; new_row];
        writetable(Rslts, fullfile(Fitdir, 'BestRslts.txt'), 'Delimiter', '\t');
        fprintf('Subject %d, model %i, nll = %f\n', subj, modeli, fval);
    end
end
delete(gcp("nocreate"));
%% 
function nll = McFadden(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); %exp(x(2))-1; % late noise difference between time-pressure conditions
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
tmp = squeeze(sum(max_from_each_distribution, 1));
probs = tmp ./ sum(tmp,2);
nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
if gpuparallel
    nll = gather(nll);
end
end

function nll = Mdl2(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); %exp(x(2))-1; % late noise difference between time-pressure conditions
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
tmp = squeeze(sum(max_from_each_distribution, 1));
probs = tmp ./ sum(tmp,2);
nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
if gpuparallel
    nll = gather(nll);
end
end

function nll = DNM(x, dat)
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
Mp = x(1); % change M to be M', absorbing the magnitude of late noise
eta = 1; % after the transformation, the late noise term is standardized as 1
delta = x(2); %exp(x(2))-1; % late noise difference between time-pressure conditions
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
tmp = squeeze(sum(max_from_each_distribution, 1));
probs = tmp ./ sum(tmp,2);
nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
if gpuparallel
    nll = gather(nll);
end
end

function nll = dnDNM(x, dat, mode) % cut inputs, independent
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
delta = x(2); %exp(x(2))-1; % late noise difference between time-pressure conditions
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
tmp = squeeze(sum(max_from_each_distribution, 1));
probs = tmp ./ sum(tmp,2);
nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
if gpuparallel
    nll = gather(nll);
end
end

function nll = dnDNMe(x, dat, mode) % cut inputs, independent
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
delta = x(2); %exp(x(2))-1; % late noise difference between time-pressure conditions
wp = x(3); % do the same transformation on w
scl = x(4); % scaling parameter on the early noise
delta2 = x(5); %exp(x(5))-1; % possible time pressure impact on the scaling parameter of early noise
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3','chosenItem','TimeConstraint'});
num_samples = 20000;
Ntrl = size(dat,1);
if strcmp(mode, 'absorb')
    samples = [];
    for ci = 1:3
        if gpuparallel
            values = gpuArray(data.(['V',num2str(ci)])');
            stds = gpuArray(data.(['sdV', num2str(ci)])').*(1 + delta2*(data.TimeConstraint'==1.5))*scl;
            samples(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        else
            values = data.(['V',num2str(ci)])';
            stds = data.(['sdV', num2str(ci)])'.*(1 + delta2*(data.TimeConstraint'==1.5))*scl;
            samples(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        end
    end
    D1 = [];
    D2 = [];
    D3 = [];
    for ci = 1:3
        if gpuparallel
            values = gpuArray(data.(['V',num2str(ci)])');
            stds = gpuArray(data.(['sdV', num2str(ci)])').*(1 + delta2*(data.TimeConstraint'==1.5))*scl;
            D1(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D2(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            D3(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        else
            values = data.(['V',num2str(ci)])';
            stds = data.(['sdV', num2str(ci)])'.*(1 + delta2*(data.TimeConstraint'==1.5))*scl;
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
tmp = squeeze(sum(max_from_each_distribution, 1));
probs = tmp ./ sum(tmp,2);
nll = -sum(log(max(probs(sub2ind(size(probs), 1:size(probs, 1), choice)), eps)));
if gpuparallel
    nll = gather(nll);
end
end
