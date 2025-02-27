function [DsnMtx]=Sampling(subNo)
%% searching bid record
logtxt = fullfile('log','txtDat',strcat('BidTask_',num2str(subNo),'.txt'));
if ~exist(logtxt,'file')
    error('Caution!! subNo %i bid record not found',subNo);
else
    bid = tdfread(logtxt);
end

%% loading BDM auction data
items = unique(bid.item);
auction = grpstats(bid.bid,{bid.item},'mean');
patch = grpstats(bid.patch,{bid.item},'mean');
group = mean(bid.Group);
%% handle output
outdir = fullfile('log','Sampling');
plotdir = fullfile(outdir,'QualityCheck');
if ~exist(outdir,'dir')
    mkdir(outdir);
    mkdir(plotdir);
end
DatPckg_dir = fullfile('DataPackages',['DatPckg_', num2str(subNo)]);
if ~exist(DatPckg_dir,'dir')
    mkdir(DatPckg_dir);
end
fontsize = 14;
mksz = 8;
aspect = [9, 3.8];
aspect1 = [9,7];
aspect2 = [4.5, 3.8];
aspect3 = [9,3.8];

%% auction data  quality check
% Pre-quality check on Auction value, give warning when
% a. warining when bidding variance is too small
if var(auction) < 100
    warning('Auction data Variance too small, please check distribution');
end

% Item sampling
% For each patch:
% 1. to take the higest NGrid number of items from the sorting as targets
% (precise patch only, and discard the vague patch high-value items)
% 2. to pick 1,3,5,7,9,11 and 2,4,6,8,10,12 ranked item for distractors,
% randomly assigned as high and low time pressure conditions.

[sauction, I] = sort(auction);
sitems = items(I);
spatch = patch(I);

NGrid = 6;
h = figure;
filename = sprintf('Bid_and_Pick_%i', subNo);
subplot(1,2,1); hold on;
H = histogram(auction, 3*NGrid);
xlabel('BDM Auction ($)');
ylabel('Frequency');
title('All items');
savefigs(h, filename, plotdir, fontsize, aspect);

subplot(1,2,2); hold on;
prsval = sauction(spatch == group);
lprs = plot(prsval,prsval,'r.','MarkerSize',mksz);
plot(prsval(end-NGrid+1:end),prsval(end-NGrid+1:end),'ko','MarkerSize',mksz);
vgval = sauction(spatch ~= group);
lvg = plot(vgval(1:(end-NGrid)),vgval(1:(end-NGrid)),'b.','MarkerSize',mksz);

xlabel('BDM Auction ($)');
ylabel('BDM Auction ($)');
legend([lprs,lvg], {'Precise','Vague'},'Location','SouthEast');
savefigs(h, filename, plotdir, fontsize, aspect);
copyfile(fullfile(plotdir,[filename, '.png']),DatPckg_dir);

%% generating choice trials
% combination of NGrid high precise items
prsitem = sitems(spatch == group);
vgitem = sitems(spatch ~= group);
C = nchoosek(prsitem(end-NGrid+1:end),2);
% 4 conditions
% Low Time pressure, decision noise coded as 1
% - precise, representation noise coded as 1
indxLP = randi(2);
interlevedP = [1:2:(numel(prsitem)-NGrid); 2:2:(numel(prsitem)-NGrid)]';
lowitem = prsitem(interlevedP(:,indxLP));
higharray = repmat(C, numel(lowitem),1);
lowarray = repmat(lowitem,1,length(C))';
IDs = [higharray, lowarray(:)];
locIDs = Shuffle(IDs')';
%           V1, V2, V3  VL, VD, VR,    decision noise code     representation noise code    generating order
LPcombo = [IDs, locIDs, repmat(1,numel(lowarray),1), repmat(1,numel(lowarray),1), [1:numel(lowarray)]'];
% - vague, representation noise coded as 2
indxLV = randi(2);
interlevedV = [1:2:(numel(vgitem)-NGrid); 2:2:(numel(vgitem)-NGrid)]';
lowitem = vgitem(interlevedV(:,indxLV));
higharray = repmat(C, numel(lowitem),1);
lowarray = repmat(lowitem,1,length(C))';
IDs = [higharray, lowarray(:)];
locIDs = Shuffle(IDs')';
%           V1, V2, V3  VL, VD, VR,    decision noise code      representation noise code    generating order
LVcombo = [IDs, locIDs, repmat(1,numel(lowarray),1), repmat(2,numel(lowarray),1), numel(lowarray)+[1:numel(lowarray)]'];
% High Time pressure, decision noise code as 2
% - precise, representation noise code as 1
lowitem = prsitem(interlevedP(:,3-indxLP));
higharray = repmat(C, numel(lowitem),1);
lowarray = repmat(lowitem,1,length(C))';
IDs = [higharray, lowarray(:)];
locIDs = Shuffle(IDs')';
%           V1, V2, V3  VL, VD, VR,    decision noise code      representation noise code    generating order
HPcombo = [IDs, locIDs, repmat(2,numel(lowarray),1), repmat(1,numel(lowarray),1), 2*numel(lowarray)+[1:numel(lowarray)]'];
% - vague, representation noise code as 2
lowitem = vgitem(interlevedV(:,3-indxLV));
higharray = repmat(C, numel(lowitem),1);
lowarray = repmat(lowitem,1,length(C))';
IDs = [higharray, lowarray(:)];
locIDs = Shuffle(IDs')';
%           V1, V2, V3  VL, VD, VR,    decision noise code      representation noise code    generating order
HVcombo = [IDs, locIDs, repmat(2,numel(lowarray),1), repmat(2,numel(lowarray),1), 3*numel(lowarray)+[1:numel(lowarray)]'];
Lcombo = Shuffle([LPcombo;LVcombo],2);
Lcombo = Shuffle(Lcombo,2); % shuffle twice to generate more randomization
Hcombo = Shuffle([HPcombo;HVcombo],2);
Hcombo = Shuffle(Hcombo,2);
ord = 2-mod(subNo,2);
switch ord
    case 1
        combo = [Lcombo; Hcombo];
    case 2
        combo = [Hcombo; Lcombo];
end

DsnMtx = table(combo(:,1),combo(:,2),combo(:,3),...
    auction(combo(:,1)),auction(combo(:,2)),auction(combo(:,3)),...
    combo(:,4),combo(:,5),combo(:,6),...
    auction(combo(:,4)),auction(combo(:,5)),auction(combo(:,6)),...
    combo(:,7),combo(:,8),combo(:,9),[1:length(combo)]',...
    'VariableNames',{'ID1','ID2','ID3','bd1','bd2','bd3','IDL','IDD','IDR','bdL','bdD','bdR','TimePressure','Vagueness','GenOrder','Trial'});
filename = sprintf('DsnMtx_%i',subNo);
save(fullfile(outdir,[filename '.mat']),'DsnMtx');
copyfile(fullfile(outdir,[filename '.mat']),DatPckg_dir);
%% Quality checks
% Spatial
h = figure;
filename = sprintf('DesignCheck_%i', subNo);
subplot(2,2,1); hold on;
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 1;
lLP = plot3(DsnMtx.bd1(mask),DsnMtx.bd2(mask),DsnMtx.bd3(mask),'c.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 2;
lLV = plot3(DsnMtx.bd1(mask),DsnMtx.bd2(mask),DsnMtx.bd3(mask),'b.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 1;
lHP = plot3(DsnMtx.bd1(mask),DsnMtx.bd2(mask),DsnMtx.bd3(mask),'m.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 2;
lHV = plot3(DsnMtx.bd1(mask),DsnMtx.bd2(mask),DsnMtx.bd3(mask),'r.','MarkerSize',mksz);
xlabel('V1');
ylabel('V2');
zlabel('V3');
grid on;
view([30,30]);
savefigs(h, filename, plotdir, fontsize, aspect1);

subplot(2,2,2); hold on;
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 1;
plot3(DsnMtx.bdL(mask),DsnMtx.bdD(mask),DsnMtx.bdR(mask),'c.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 2;
plot3(DsnMtx.bdL(mask),DsnMtx.bdD(mask),DsnMtx.bdR(mask),'b.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 1;
plot3(DsnMtx.bdL(mask),DsnMtx.bdD(mask),DsnMtx.bdR(mask),'m.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 2;
plot3(DsnMtx.bdL(mask),DsnMtx.bdD(mask),DsnMtx.bdR(mask),'r.','MarkerSize',mksz);
xlabel('VL');
ylabel('VD');
zlabel('VR');
grid on;
view([30,30]);
savefigs(h, filename, plotdir, fontsize, aspect1);
% Temporal
subplot(2,1,2); hold on;
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 1;
lLP = plot(DsnMtx.Trial(mask),DsnMtx.bd1(mask),'c.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd2(mask),'c.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd3(mask),'c.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 1 & DsnMtx.Vagueness == 2;
lLV = plot(DsnMtx.Trial(mask),DsnMtx.bd1(mask),'b.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd2(mask),'b.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd3(mask),'b.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 1;
lHP = plot(DsnMtx.Trial(mask),DsnMtx.bd1(mask),'m.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd2(mask),'m.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd3(mask),'m.','MarkerSize',mksz);
mask = DsnMtx.TimePressure == 2 & DsnMtx.Vagueness == 2;
lHV = plot(DsnMtx.Trial(mask),DsnMtx.bd1(mask),'r.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd2(mask),'r.','MarkerSize',mksz);
plot(DsnMtx.Trial(mask),DsnMtx.bd3(mask),'r.','MarkerSize',mksz);
legend([lLP,lLV,lHP,lHV],{'Low Precise','Low Vague','High Precise','High Vague'},'FontSize',fontsize-5);
xlabel('Trial');
ylabel('Bid');
savefigs(h, filename, plotdir, fontsize, aspect1);
end
