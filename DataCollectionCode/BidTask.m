%function BidTask(subNo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The bid task outside of scanner
% Created by Bo Shen, NYU School of Medicine, Jun/14/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------- Input Dialogue ---------------
if 1%nargin < 1
    % default settings, subID string [Year Month Day Hour Min Sec]
    if ~exist('subNo','var')
        subID = datestr(now,30);
        subID([1,2,9,14,15])=[];
        group = '1';
        promptParameters = {'Subject ID', 'Group'};
        defaultParameters = {subID,group};
        Settings = inputdlg(promptParameters, 'Settings', 1,  defaultParameters);
        subID = Settings{1};
        subNo = str2double(subID);
        group = str2double(Settings{2});
    end
end
addpath('func');
%% ------- Initialize MATLAB ---------------
rand('state',sum(100*clock)); % reset random seed
KbName('UnifyKeyNames'); % Unify Keyboard Across Platform
[~, ~, KeyCode] = KbCheck; % initializing functions
KeyCode(KeyCode==1) = 0;
GetSecs;

%% ------- Output file handling ---------------
log_dir = fullfile('log');
% creat log folders if not exist
if ~exist(log_dir,'dir')
    mkdir(log_dir);
    mkdir(fullfile(log_dir,'txtDat'));
end
DatPckg_dir = fullfile('DataPackages',['DatPckg_', num2str(subNo)]);
if ~exist(DatPckg_dir,'dir')
    mkdir(DatPckg_dir);
end
scrn_dir = fullfile(log_dir,'Screenshot',num2str(subNo));
% if ~exist(scrn_dir,'dir')
%     mkdir(scrn_dir);
% end
logtxt = fullfile(log_dir,'txtDat',strcat('BidTask_',num2str(subNo),'.txt'));
if fopen(logtxt, 'rt')~=-1
    fclose('all');
    error('Caution! This Subject ID already exists! Choose a different Subject ID.');
else
    datafilepointer = fopen(logtxt,'wt'); % open ASCII file for writing
end
fprintf(datafilepointer,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Group','subID','trial','item','itemname','patch','bid','bid_times','certainty',...
    'fixation_time','bid_start_time','bid_finish_time','RT');
try
    %% ------- Setting PTB Display ---------------
    AssertOpenGL;
    screens = Screen('Screens');
    if length(screens) == 1
        ScrnNo = 0;
    elseif length(screens) >= 2 % in just having external monitors, choose the main
        ScrnNo = 1;% min(screens);
    end
    backgroundgrayindex = 255/255;  % uncomment to set background lumiNaNce 0 - 1
    if exist('backgroundgrayindex','var') % in case manually set gray index
        background = GrayIndex(ScrnNo,backgroundgrayindex);
    else
        background = GrayIndex(ScrnNo);
    end
    Screen('Preference', 'SkipSyncTests', 1);
    [Win, wRect] = Screen('OpenWindow',ScrnNo, background);
    [xcenter, ycenter] = RectCenter(wRect);
    fps = Screen('FrameRate',Win); % frames per second
    ifi = Screen('GetFlipInterval', Win); % inter flip interval
    if fps == 0
        fps = 1/ifi;
    end
    HideCursor;	% Hide the mouse cursor
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', Win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    % Set priority for script execution to realtime priority:
    priorityLevel=MaxPriority(Win);
    Priority(priorityLevel);
    
    %% ------- Setting for Display Contents ---------------
    % Colorspace
    blue = [0,153,255];
    orange = [255,125,0];
    green = [0,255,0];
    red = [255,0,0];
    white = [255,255,255];
    black = [0 0 0];
    gray = [127 127 127];
    yellow = [255 255 0];
    % font control
    Screen('TextFont',Win, 'Arial');
    font_size = 28;
    smallfont_size = 24;
    % Display positions and frames
    itemPosition = [xcenter, ycenter, xcenter, ycenter];
    namebox = [xcenter-400, itemPosition(4)+font_size/2, xcenter+400, itemPosition(4)+3/2*font_size];
    widst = round(wRect(3)*.8);%1000;
    scalebox = [xcenter-widst/2, namebox(4)+3*font_size, xcenter+widst/2, namebox(4)+4*font_size];
    VASycenter = mean(scalebox([2,4]));
    VASHeight = scalebox(4) - scalebox(2) + 2*font_size;
    scaletext = 'The highest $ willing to pay for it?';
    scaleTitle = [scalebox(1), scalebox(2)-font_size*3/2, scalebox(3), scalebox(2)-font_size/2];
    scale = {'$0','$90'};
    scaleboxL = [scalebox(1), scalebox(4)+font_size/2, scalebox(1)+font_size, scalebox(4)+3/2*font_size];
    scaleboxR = [scalebox(3)-font_size, scalebox(4)+font_size/2, scalebox(3), scalebox(4)+3/2*font_size];
    
    certaintybox = [xcenter-round(widst/3), scaleboxL(4)+font_size, xcenter+round(widst/3), scaleboxL(4)+2*font_size];
    VAS2ycenter = mean(certaintybox([2,4]));
    VAS2Height = certaintybox(4) - certaintybox(2) + 2*font_size;
    certaintext = 'How certain are you?';
    certainTitle = [certaintybox(1), certaintybox(2)-font_size*3/2, certaintybox(3), certaintybox(2)-font_size/2];
    certainboxL = [certaintybox(1), certaintybox(4)+font_size/2, certaintybox(1)+font_size, certaintybox(4)+3/2*font_size];
    certainboxR = [certaintybox(3)-font_size, certaintybox(4)+font_size/2, certaintybox(3), certaintybox(4)+3/2*font_size];
    certain = {'Extremely uncertain: 0','Absolutely certain: 100'};
    checkbox = [xcenter+2*font_size, scaleboxL(4)+5/4*font_size, xcenter+4*font_size, scaleboxL(4)+9/4*font_size];
    Nextbox = [xcenter - 3/2*font_size, checkbox(4) + 2*font_size, xcenter + 3/2*font_size, checkbox(4) + 4*font_size];
    
    %% ------- Define response bottons --------------
    enter_button = KbName('Return'); % controled by experimenter
    escapeKey = KbName('`~'); % escape button, controlled by experimenter
    [x, y, mouse] = GetMouse(Win); % sliding and confirm button
    
    %% ------- Import Images and icons ---------------
    [num, txt, raw] = xlsread(fullfile('.','CorrectStimuli','CorrectItems.xlsx'), 1, 'A:F');
    Nimgs = length(num);
    itemNames = txt(2:end,2);
    [Y, I] = sort(num(:,1));
    num = num(I,:);
    itemNames = itemNames(I);
    
    itemlist = dir(fullfile('.','CorrectStimuli','*.*g'));
    itemlist(strcmp({itemlist.name},'Fixation.jpg')) = [];
    for i = 1:length(itemlist)
        namelist{i} = itemlist(i).name(1:end-4);
    end
    
    check = [];
    for i = 1:length(itemNames)
        check(i) = any(strcmpi(itemNames{i},namelist));
    end
    if any(~check)
        error('The name unmatached between the Spreadsheet and the Stimuli folder');
    end
    
    for item = 1:Nimgs % Full items
        name = itemNames{item};
        myimgfile = dir(fullfile('.','CorrectStimuli',sprintf('%s.*',lower(name))));
        img = imread(fullfile(myimgfile.folder, myimgfile.name));
        imgaspect(item,:) = [size(img,2),size(img,1)];
        itemTx{item} = Screen('MakeTexture', Win, img);
    end
    DiagLen = sqrt(imgaspect(:,1).^2 + imgaspect(:,2).^2);
    PicInch = 400;
    ImgInch = [imgaspect(:,1)./DiagLen imgaspect(:,2)./DiagLen]*PicInch; % normalize the image to be the same diagnal inches
    ImgBox = round([-ImgInch(:,1)/2, -ImgInch(:,2), ImgInch(:,1)/2, 0*ImgInch(:,2)]);
    Ntrials = 3*Nimgs;
    %% ------- Experiment Begin! ---------------
    % Introduction
    HideCursor;
    Screen('TextSize', Win, font_size);
    Instruction = ['Please bid for each item you see in the next rounds of the task\n'...
                    'Remember, how much you get depends on how good your decisions are\n'...
                    'Press any key to continue'];
    DrawFormattedText(Win,double(Instruction),'center','center', black,[],[],[],2);
    Screen('Flip', Win);
    % Waiting for any key press
    t.enter = KbStrokeWait;
    t.bidfinish = t.enter;
    % Bidding start...
    list = [Shuffle(1:Nimgs) Shuffle(1:Nimgs) Shuffle(1:Nimgs)]; % randomize the order
    while any(list(1:end-1) == list(2:end)) % make sure no adjacent same item
        list = [Shuffle(1:Nimgs) Shuffle(1:Nimgs) Shuffle(1:Nimgs)];
    end
    for trial = 1:Ntrials
        bid_times = ceil(trial/Nimgs);
        item = list(trial);
        patch = num(num(:,1) == item,6);
        DrawFormattedText(Win,double('+'),'center','center', black,[],[],[],2);
        t.fix = Screen('Flip', Win, t.bidfinish + .5);
        ITI =  .5+rand(1,1)*2.5; % uniform distrib. from .5 to 3 secs
        name = itemNames{item};
        t.bidstart = Screen('Flip', Win, t.fix + ITI);
        % SetMouse(cx,cy,Win);
        ShowCursor('Arrow');
        bid = NaN;
        rate = NaN;
        next = 0;
        tic;
        while isnan(bid) || next == 0
            [x, y, mouse] = GetMouse(Win);
            [~, ~, KeyCode] = KbCheck;
            Escape(KeyCode, escapeKey);
            
            Screen('DrawTexture', Win, itemTx{item}, [], itemPosition+ImgBox(item,:)); % item picture
            DrawFormattedText(Win,double(name),'center','center',black,[],[],[],1,[],namebox); % item name
            DrawFormattedText(Win,double(scaletext),'center','center',black,[],[],[],1,[],scaleTitle); % bidding instruction
            Screen('FillRect', Win, gray, scalebox); % slider bar
            DrawFormattedText(Win,double(scale{1}),'center','center',black,[],[],[],1,[],scaleboxL); % price ticker at left, $0
            DrawFormattedText(Win,double(scale{2}),'center','center',black,[],[],[],1,[],scaleboxR); % price ticker at right, $90
            DrawFormattedText(Win,double(certaintext),'center','center',black,[],[],[],1,[],certainTitle); % certainty rating instruction
            Screen('FillRect', Win, gray, certaintybox); % uncertainty rating bar
            DrawFormattedText(Win,double(certain{1}),'center','center',black,[],[],[],1,[],certainboxL); % rating to left, absolutely sure
            DrawFormattedText(Win,double(certain{2}),'center','center',black,[],[],[],1,[],certainboxR); % extremely uncertain
            
            Screen('FillRect', Win, yellow, Nextbox); % next button
            DrawFormattedText(Win,double('>>'),'center','center',black,[],[],[],1.4,[],Nextbox+[0,-font_size/4,0,-font_size/4]);
            
            if (y >= VASycenter-VASHeight/2 && y <= VASycenter+VASHeight/2) % IsInRect(x,y,scalebox)
                if (x < xcenter - widst/2)
                    x = xcenter - widst/2;
                elseif (x > xcenter + widst/2)
                    x = xcenter + widst/2;
                end
                ShowCursor('CrossHair');
                Screen('DrawLines', Win, [x x; VASycenter-font_size VASycenter+font_size], 1, black,[0 0],1);
                value = round((x - xcenter + widst/2)/widst*90,1);
                valuebox = [x-font_size VASycenter-font_size x+font_size VASycenter+font_size/1.4];
                DrawFormattedText(Win,double(sprintf('%.1f',value)),'center','center',black,[],[],[],1.4,[],valuebox);
                if any(mouse)
                    bid = value;
                    bidbox = valuebox;
                end
            elseif (y >= VAS2ycenter-VAS2Height/2 && y <= VAS2ycenter+VAS2Height/2) && ~isnan(bid) % IsInRect(x,y,checkbox)
                if (x < xcenter - widst/3)
                    x = xcenter - widst/3;
                elseif (x > xcenter + widst/3)
                    x = xcenter + widst/3;
                end
                ShowCursor('CrossHair');
                Screen('DrawLines', Win, [x x; VAS2ycenter-font_size VAS2ycenter+font_size], 1, black,[0 0],1);
                value = round((x - xcenter + widst/3)/widst/2*3*100);
                valuebox = [x-font_size VAS2ycenter-font_size x+font_size VAS2ycenter+font_size/1.4];
                DrawFormattedText(Win,double(sprintf('%i',value)),'center','center',black,[],[],[],1.4,[],valuebox);
                if any(mouse)
                    rate = value;
                    ratebox = valuebox;
                end
            elseif IsInRect(x,y,Nextbox) && ~isnan(bid)
                ShowCursor('Hand');
                if any(mouse)
                    next = 1;
                end
            else
                ShowCursor('Arrow');
            end
            if ~isnan(bid)
                DrawFormattedText(Win,double(sprintf('%.1f',bid)),'center','center',red,[],[],[],1.4,[],bidbox);
            end
            if ~isnan(rate)
                DrawFormattedText(Win,double(sprintf('%i',rate)),'center','center',red,[],[],[],1.4,[],ratebox);
            end
            t.bidfinish = Screen('Flip', Win);
        end
        %imgfile =  fullfile(scrn_dir,sprintf('BidTrial_%03i.jpg',trial));
        %ImgArray = Screen('GetImage',Win);
        %imwrite(ImgArray,imgfile);
        %WaitSecs(1.0);
        fprintf(datafilepointer,'%i\t%s\t%i\t%i\t%s\t%i\t%.1f\t%i\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            group,subID,trial,item,name,patch,bid,bid_times,rate,...
            t.fix-t.enter,t.bidstart-t.enter,t.bidfinish-t.enter,t.bidfinish-t.bidstart);
        copyfile(logtxt,DatPckg_dir);
    end
    message = 'The current task is complete.\nPlease wait';
    DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
    Screen('Flip', Win,[],1);
catch
    % catch error: This is executed in case something goes wrong in the
    % 'try' part due to programming error etc.:
    % Do same cleanup as at the end of a regular session...
    sca;
    ShowCursor;
    fclose('all');
    Priority(0);
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end

%% Processing
message = '\n\n\n Sampling...';
DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
Screen('Flip', Win,[],1);
Sampling(subNo);


% %% Data transfer
% message = '\n\n\n\n\n\n\n Transferring...';
% DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
% Screen('Flip', Win);
% Indvfolder = sprintf('DataTransfer_%i',subNo);
% mkdir(Indvfolder);
% copyfile(logtxt,Indvfolder);
% DsnMtx = fullfile('Sampling',sprintf('DsnMtx_%i.mat',subNo));
% copyfile(DsnMtx,Indvfolder);
% Seq = fullfile('Sequence',sprintf('Sequence_%i.mat',subNo));
% copyfile(Seq,Indvfolder);

sca;
ShowCursor;
fclose('all');
Priority(0);
%% Escape
function Escape(KeyCode, escapeKey)
if ( KeyCode(escapeKey) == 1 )
    Priority(0);
    sca;
    ShowCursor;
    fclose('all');
    error('experiment aborted by user');
end
end
