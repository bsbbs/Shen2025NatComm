%% Load in raw data and do preliminary process
% Define I/O directories
datadir = fullfile('.', 'myData','RawData');
svdir = fullfile('.', 'myData');

bd = load_data(datadir, 'BidTask_22*');
bd.Definitive = (bd.Group == bd.patch);
bd.Definitive = double(bd.Definitive);
% Rename the variables
bd = renamevars(bd, {'item'}, {'Item'});
disp(head(bd));

%%
% Raw scale of bid value
bdw = unstack(bd(:,[1,2,4,6,7,8,14]), 'bid', 'bid_times');
bdw.BidMean = mean(bdw(:, ["x1", "x2", "x3"]).Variables, 2);
bdw.BidSd = std(bdw(:, ["x1", "x2", "x3"]).Variables, 0, 2);
bdw.sd12 = std(bdw(:, ["x1", "x2"]).Variables, 0, 2);
bdw.sd23 = std(bdw(:, ["x2", "x3"]).Variables, 0, 2);
bdw.sd13 = std(bdw(:, ["x1", "x3"]).Variables, 0, 2);
bdw.cv = bdw.BidSd ./ bdw.BidMean;
disp(head(bdw));
%% To load choice data
mt = load_data(datadir, 'MainTask_22*');
mt.Vaguenesscode = mt.Vaguenesscode - 1;
mt.Definitive = 1 - mt.Vaguenesscode;
disp(head(mt));
%% To merge the bidding variance information into the choice matrix
IDs = {'ID1', 'V1'; 'ID2', 'V2'; 'ID3', 'V3'};
for i = 1:size(IDs, 1)
    ID = IDs{i, 1};
    V = IDs{i, 2};
    mt.Item = mt.(ID);
    tmp = innerjoin(bdw(:, {'subID', 'Item', 'BidSd'}), mt(:, {'subID', 'Item', 'trial'}));
    tmpp = renamevars(tmp, 'BidSd', ['sd' V]);
    mt = innerjoin(mt, tmpp(:, {'subID', 'trial', ['sd' V]}));
end
mt = removevars(mt, 'Item');
mt = rmmissing(mt, 'DataVariables', 'chosenItem');
mt.chosenItem = int32(mt.chosenItem);
disp(head(mt));
%% save the preprocessed data
save(fullfile(svdir, 'TrnsfrmData.mat'), "mt");
writetable(mt, fullfile(svdir, 'TrnsfrmData.csv'));
%%
function df = load_data(directory, file_pattern)
file_paths = dir(fullfile(directory, file_pattern));
data = cell(numel(file_paths), 1);

for i = 1:numel(file_paths)
    file_path = fullfile(file_paths(i).folder, file_paths(i).name);
    file_data = readtable(file_path, 'Delimiter', '\t');
    data{i} = file_data;
end

df = vertcat(data{:});
end
