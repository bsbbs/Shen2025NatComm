function [probs, Ovlps, CVs] = dnDNM(dat, pars, constraint, products)
%%%%%%%%%%%%%%%%%%%%%%%
% Dual-noise Divisive Normalization Model
% Created by Bo Shen, NYU, 2024
%
% Input structure:
% - dat: data structured as mean values of V1, V2, V3, and their standard deviation sdV1, sdV2, and sdV3
% - pars: controlling parameters, eta, scl, w, and M
% - constraint: implementing biological/non-negative constraint or not
% - products: the requested variables of calculation for outputs
%%%%%%%%%%%%%%%%%%%%%%%%
%% detect if you have NVIDA GPU device to accellarate
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
%% parameters to define
num_samples = 1024*1e3; % number of replication for each set
data = dat(:, {'V1', 'V2', 'V3', 'sdV1','sdV2','sdV3'}); % data comes in with mean values and their early noise std
Ntrl = size(dat,1); % number of the testing trials
eta1 = pars(1); % late noise standard deviation
eta2 = pars(2); % late noise standard deviation
eta3 = pars(3); % late noise standard deviation
scl = pars(4); % scaling for early noise
w = pars(5); % weight of normalization
M = pars(6); % baseline normalization
K = 75; % the neural maximum response as assumed
%% simulation begin...
Numerator = [];
Denominator1 = [];
Denominator2 = [];
Denominator3 = [];
for ci = 1:3
    if gpuparallel
        values = gpuArray(data.(['V',num2str(ci)])');
        stds = gpuArray(data.(['sdV', num2str(ci)])')*scl;
        if strcmp(constraint, 'none')
            Numerator(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator1(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator2(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator3(:,:,ci) = gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        elseif strcmp(constraint, 'biological')
            Numerator(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator1(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator2(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator3(:,:,ci) = max(gpuArray.randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        end
    else
        values = data.(['V',num2str(ci)])';
        stds = data.(['sdV', num2str(ci)])'*scl;
        if strcmp(constraint, 'none')
            Numerator(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator1(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator2(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
            Denominator3(:,:,ci) = randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1);
        elseif strcmp(constraint, 'biological')
            Numerator(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator1(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator2(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
            Denominator3(:,:,ci) = max(randn([num_samples, Ntrl]).*stds + repmat(values, num_samples, 1), 0);
        end
    end
end
Denominator1 = sum(Denominator1, 3)*w + M;
Denominator2 = sum(Denominator2, 3)*w + M;
Denominator3 = sum(Denominator3, 3)*w + M;
% The product of divisive normalization before adding late noise
DNP = K*Numerator./cat(3, Denominator1, Denominator2, Denominator3);
clear Denominator1 Denominator2 Denominator3;
if gpuparallel
    if strcmp(constraint, 'none')
        SVs = DNP + cat(3, gpuArray.randn([num_samples, Ntrl])*eta1, gpuArray.randn([num_samples, Ntrl])*eta2, gpuArray.randn([num_samples, Ntrl])*eta3);
    elseif strcmp(constraint, 'biological')
        SVs = max(DNP + cat(3, gpuArray.randn([num_samples, Ntrl])*eta1, gpuArray.randn([num_samples, Ntrl])*eta2, gpuArray.randn([num_samples, Ntrl])*eta3), 0);
    end
else
    if strcmp(constraint, 'none')
        SVs = DNP + cat(3, randn([num_samples, Ntrl])*eta1, randn([num_samples, Ntrl])*eta2, randn([num_samples, Ntrl])*eta3);
    elseif strcmp(constraint, 'biological')
        SVs = max(DNP + cat(3, randn([num_samples, Ntrl])*eta1, randn([num_samples, Ntrl])*eta2, randn([num_samples, Ntrl])*eta3), 0);
    end
end
clear DNP;
%% calculating outputs
if ismember('Probability', products)
    probs = CalculateProbs(SVs);
    if gpuparallel
        probs = gather(probs);
    end
else
    probs = NaN;
end
if ismember('Coeff_of_Var', products)
    CVs = squeeze(std(SVs, [], 1)./mean(SVs, 1));
    if gpuparallel
        CVs = gather(CVs);
    end
else
    CVs = NaN;
end
if ismember('Overlap', products)
    Ovlps = nan([1, numel(dat.V3)]);
    if gpuparallel
        SVs = gather(SVs);
    end
    dSVrng = [min(SVs(:)), max(SVs(:))];
    for v3i = 1:numel(dat.V3)
        pd1 = fitdist(SVs(:,v3i, 1),'kernel','Kernel','normal');
        x = dSVrng(1):.1:dSVrng(2);
        y1 = pdf(pd1, x);
        pd2 = fitdist(SVs(:,v3i, 2),'kernel','Kernel','normal');
        y2 = pdf(pd2, x);
        ovlp = sum([y2(y1>y2), y1(y2>y1)])*.1;
        Ovlps(v3i) = ovlp*100;
    end
else
    Ovlps = NaN;
end
end