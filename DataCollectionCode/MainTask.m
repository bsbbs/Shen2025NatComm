function MainTask(subNo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main task in fMRI scanning
% Created by Bo Shen, NYU School of Medicine, May/5/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------- Input Dialogue ---------------
if nargin < 1
    % default settings, subID string [Year Month Day Hour Min Sec]
    if ~exist('subNo','var')
        subID = datestr(now,30);
        subID(9)=[];
        promptParameters = {'Subject ID'};
        defaultParameters = {subID};
        Settings = inputdlg(promptParameters, 'Settings', 1,  defaultParameters);
        subID = Settings{1};
        subNo = str2double(subID);
    end
end
subID = num2str(subNo);
%% ------- Loading Design Matrix and Sequence
DatPckg_dir = fullfile('DataPackages',['DatPckg_', num2str(subNo)]);
try
    load(fullfile(DatPckg_dir,sprintf('DsnMtx_%s.mat',num2str(subNo))));
catch
    error('Sampling file is missing for Subject ID %i', subNo);
end

%% ------- Define experiment parameters ----
Nptrials = 14; % number of practicing trials
including_practice = 1;
Ntrial = length(DsnMtx.Trial);
breakat = 1:60:(Ntrial+1);
breakat(end) = [];
TPvalue = [10, 1.5];
TPcode = {'Low','High'};
Vaguecode = {'Precise','Vague'};
%% ------- Initialize MATLAB ---------------
rand('state',sum(100*clock)); % reset random seed
KbName('UnifyKeyNames'); % Unify Keyboard Across Platform
[~, ~, KeyCode] = KbCheck; % initializing functions
KeyCode(KeyCode==1) = 0;
GetSecs;

%% ------- Output file handling ---------------
log_dir = fullfile('log');
txt_dir = fullfile(log_dir,'txtDat');
scrn_dir = fullfile(log_dir,'Screenshot',num2str(subNo));
% creat log folders if not exist
if ~exist(log_dir,'dir')
    mkdir(log_dir);
    mkdir(txt_dir);
    mkdir(fullfile(log_dir,'Screenshot'));
end
% if ~exist(scrn_dir,'dir')
%     mkdir(scrn_dir);
% end
logtxt = fullfile(txt_dir,strcat('MainTask_',num2str(subNo),'.txt')); % logfile for the main task
practxt = fullfile(txt_dir,strcat('Practice_',num2str(subNo),'.txt')); % logfile for the practice task
if fopen(logtxt, 'rt')~=-1
    fclose('all');
    error('Caution! This Subject ID already exists! Choose a different Subject ID.');
else
    datafilepointer = fopen(logtxt,'wt'); % open ASCII file for writing
    datafilepointerp = fopen(practxt,'wt');
end

fprintf(datafilepointer,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
            'Group','subID','trial', 'TimePressure', 'TimeConstraint', 'Vagueness', 'Vaguenesscode',... %7
        'ID1', 'ID2', 'ID3', 'V1', 'V2', 'V3',... % 13
        'ID_left', 'ID_middle','ID_right', 'Bid_left', 'Bid_middle', 'Bid_right',... %19
        'chosenItem', 'RT', 'RTvariance', 'chosenPosition', 'chosenID', 'chosenValue',...%25
            'fixation_time', 'option_time', 'choice_time','feedback_time','timeout');%30
        
fprintf(datafilepointerp,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
                    'Group','subID','trial', 'TimePressure', 'TimeConstraint',...
                    'ID_left', 'ID_middle','ID_right',...
                    'RT', 'RTvariance', 'chosenPosition', 'chosenID',...
                    'fixation_time', 'option_time', 'choice_time','feedback_time','timeout');
try
    %% ------- Setting PTB Display ---------------
    AssertOpenGL;
    screens = Screen('Screens');
    if length(screens) == 1
        ScrnNo = 0;
    elseif length(screens) >= 2 % in just having external monitors, choose the main
        ScrnNo = max(screens);
    end
    backgroundgrayindex = 255/255;  % uncomment to set background lumiNaNce 0 - 1
    if exist('backgroundgrayindex','var') % in case manually set gray index
        background = GrayIndex(ScrnNo,backgroundgrayindex);
    else
        background = GrayIndex(ScrnNo);
    end
    Screen('Preference', 'SkipSyncTests', 1);
    [Win, wRect] = Screen('OpenWindow',ScrnNo, background); % Screen('WindowSize', ScrnNo); [1920, 1080] in scanner
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
    % Display positions and frames
    font_size = 28;
    [xcenter, ycenter] = RectCenter(wRect);
    center = [xcenter, ycenter];
    unitlen = round(wRect(4)/8);
    PicInch = unitlen;
    StandardBox = round(PicInch/sqrt(2)*[-1 -1 1 1]/2);
    ycenter = round(ycenter - unitlen/4);
    downmove = unitlen;
    uppermove = unitlen*.5;
    leftmove = unitlen*cos(pi/6);
    rightmove = unitlen*cos(pi/6);
    LeftCenter = round([xcenter-leftmove, ycenter-uppermove, xcenter-leftmove, ycenter-uppermove]);
    DownCenter = round([xcenter, ycenter+downmove, xcenter, ycenter+downmove]);
    RightCenter = round([xcenter+rightmove, ycenter-uppermove, xcenter+rightmove, ycenter-uppermove]);
    Positions = [LeftCenter;DownCenter;RightCenter]; % Position 1, 2, and 3
    FixationCenter = [xcenter, ycenter, xcenter, ycenter];
    FixationBox = FixationCenter + [-1 -1 1 1]*font_size;
    Expansion = [-1 -1 1 1]*font_size;
    linewidth = 2;
    %% ------- Define response bottons --------------
    enter_button = KbName('Return'); % controled by experimenter
    escapeKey = KbName('`~'); % escape button, controlled by experimenter
    bL = KbName('q'); % left button
    bD = KbName('space'); % down button
    bR = KbName('p'); % right button
    triggerKey = enter_button;
    buttoncode = {'L','D','R'};
    %% ------- Import Images and icons ---------------
    fixationimgfile = fullfile('.','CorrectStimuli','Fixation.jpg');
    img = imread(fixationimgfile);
    FixTx = Screen('MakeTexture', Win, img([190:319]+186*0,[1:130]+199*1,:));
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
    ImgInch = [imgaspect(:,1)./DiagLen imgaspect(:,2)./DiagLen]*PicInch; % normalize the image to be the same diagnal inches
    ImgBox = round([-ImgInch/2, ImgInch/2]);
    %% ------- Experiment Begin! ---------------
    % font control
    Screen('TextFont',Win, 'Arial');
    % Alternative Font types:'Simhei','Times','Helvetica','Geneva'
    Screen('TextSize', Win, font_size);
    if including_practice
        % Introduction
        message = 'Practice task\nPress any button to start practicing...';
        DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
        Screen('Flip', Win);
        t.practice = KbStrokeWait;
        %% Practice, using items that not picked for the main task
        %% loading BDM auction data
        bidtxt = fullfile('log','txtDat',strcat('BidTask_',num2str(subNo),'.txt'));
        if ~exist(bidtxt,'file')
            error('Caution!! subNo %i bid record not found',subNo);
        else
            bid = tdfread(bidtxt);
        end
        items = unique(bid.item);
        auction = grpstats(bid.bid,{bid.item},'mean');
        patch = grpstats(bid.patch,{bid.item},'mean');
        group = mean(bid.Group);
        [sauction, I] = sort(auction);
        sitems = items(I);
        spatch = patch(I);
        vgitem = sitems(spatch ~= group);
        NGrid = 6;
        C = nchoosek(vgitem(end-NGrid+1:end),3);
        IDs = Shuffle(C,2);
        locIDs = Shuffle(IDs')';
        
        goodtogo = 0;
        while goodtogo == 0
            TimePressure = 1;
            for ptrial = 1:Nptrials
                if ptrial == (round(Nptrials/2)+1)
                    TimePressure = 2;
                    message = sprintf('In part of the task\nYou will only have %2.1f seconds to make a decision\nMake your choice quickly!\nPress any button to experience',TPvalue(2));
                    DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
                    Screen('Flip', Win);
                    KbStrokeWait;
                end
                %% ITI Fixation
                Screen('DrawTexture', Win, FixTx, [], FixationBox);
                t.fix = Screen('Flip', Win);
                %% Choice
                IDsLDR = locIDs(ptrial,:);
                PresentBoxes = zeros(3,4);
                stacki = 0;
                for alteri = 1:3 % pre-loading stimuli
                    stacki = stacki + 1;
                    PresentBoxes(stacki,:) = Positions(alteri,:)+round(ImgBox(IDsLDR(alteri),:)*1.1);
                    Screen('DrawTexture', Win, itemTx{IDsLDR(alteri)}, [], Positions(alteri,:)+ImgBox(IDsLDR(alteri),:));
                    Screen('FrameRect',Win,black,PresentBoxes(stacki,:),linewidth);
                end
                BoxL = PresentBoxes(1,:);
                BoxD = PresentBoxes(2,:);
                BoxR = PresentBoxes(3,:);
                
                ITI = .5+rand(1,1)*2.5; % uniform distrib. from .5 to 3 secs
                t.options = Screen('Flip', Win, t.fix + ITI - ifi/2,1); % show stimuli
                
                KeyCode(KeyCode==1) = 0;
                while ~(KeyCode(bL) || KeyCode(bD) || KeyCode(bR)) && GetSecs < t.options + TPvalue(TimePressure) - ifi
                    [~, t.choice, KeyCode, RTvariance] = KbCheck;
                    Escape(KeyCode, escapeKey);
                end
                if ~(KeyCode(bL) || KeyCode(bD) || KeyCode(bR))
                    message = 'Time Out!';
                    DrawFormattedText(Win, double(message),'center','center',red,[],[],[],2);
                    t.choice = NaN;
                    RT = NaN;
                    RTvariance = NaN;
                    chosenPosition = NaN;
                    chosenID = NaN;
                    t.feedback = Screen('Flip', Win);
                    timeout = 1;
                    WaitSecs(2.0);
                else
                    choice = [KeyCode(bL) KeyCode(bD) KeyCode(bR)];
                    RT = t.choice - t.options;
                    SlctBox = choice*[BoxL; BoxD; BoxR];
                    [SlctCenterx,~] = RectCenter(SlctBox);
                    chosenPosition = buttoncode{choice==1};
                    chosenID = IDsLDR(SlctCenterx == Positions(:,1));
                    Screen('FrameRect',Win,red,SlctBox,linewidth);
                    t.feedback = Screen('Flip', Win);
                    WaitSecs(.5);
                    timeout = 0;
                end
                fprintf(datafilepointerp,'%i\t%s\t%i\t%s\t%2.1fi\t%i\t%i\t%i\t%.3f\t%.3f\t%s\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n',...
                    group, subNo, ptrial, TPcode{TimePressure}, TPvalue(TimePressure),...
                    IDsLDR(1), IDsLDR(2), IDsLDR(3),...
                    RT, RTvariance, chosenPosition, chosenID,...
                    t.fix-t.practice, t.options-t.practice, t.choice-t.practice, t.feedback-t.practice, timeout);
            end
            message = 'Ready to go?\nYes (Do the real task)                   No (Practice again)';
            DrawFormattedText(Win, double(message),'center','center',black,[],[],[],2);
            Screen('Flip', Win, t.feedback + .5 - ifi/2);
            KeyCode(KeyCode==1) = 0;
            while ~(KeyCode(bL) || KeyCode(bR))
                [~, ~, KeyCode, ~] = KbCheck;
                Escape(KeyCode, escapeKey);
            end
            if KeyCode(bL)
                goodtogo = 1;
            elseif KeyCode(bR)
                goodtogo = 0;
            end
        end
        copyfile(practxt,DatPckg_dir);
    end
    %% Main part
    % Introduction
    message = 'Remember, how much you get depends on how good your decisions are.\nPress any key to continue...';
    DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
    Screen('Flip', Win);
    WaitSecs(1.0);
    t.main = KbStrokeWait;
    
    %% enter the main task loop
    ITI = .5+rand(1,Ntrial)*2.5; % uniform distrib. from .5 to 3 secs
    TimePressure = 0;
    for trial = 1:Ntrial
        % time pressure detecting
        if TimePressure ~= DsnMtx.TimePressure(trial)
            TimePressure = DsnMtx.TimePressure(trial);
            if trial == 1
                message1 = 'In the first half of the task,';
            elseif trial == Ntrial/2+1
                message1 = 'In the rest of the task,';
            end
            if TimePressure == 1
                message2 = sprintf('\nyou will have %2.1f seconds for every decision.\nTake your time.',TPvalue(TimePressure));
            elseif TimePressure == 2
                message2 = sprintf('\nyou only have %2.1f seconds for every decision.\nMake your choice quickly.',TPvalue(TimePressure));
            end
            message3 = '\nPress any key to contibue';
            DrawFormattedText(Win,double([message1, message2, message3]),'center','center', black,[],[],[],2);
            Screen('Flip', Win);
            KbStrokeWait;
        end
        %% Session break guidance
        if any(trial == breakat)
            message = 'Here We Go!';
            DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
            Screen('Flip', Win);
            WaitSecs(2.0);
        end
        
        %% ITI Fixation
        Screen('DrawTexture', Win, FixTx, [], FixationBox);
        t.fix = Screen('Flip', Win);

        %% Choice
        IDsLDR = [DsnMtx.IDL(trial), DsnMtx.IDD(trial), DsnMtx.IDR(trial)];
        ValuesLDR = [DsnMtx.bdL(trial), DsnMtx.bdD(trial), DsnMtx.bdR(trial)];
        IDs123 = [DsnMtx.ID1(trial), DsnMtx.ID2(trial), DsnMtx.ID3(trial)];
        Values123 = [DsnMtx.bd1(trial), DsnMtx.bd2(trial), DsnMtx.bd3(trial)];
        PresentBoxes = zeros(3,4);
        stacki = 0;
        for alteri = 1:3 % pre-loading stimuli
            stacki = stacki + 1;
            PresentBoxes(stacki,:) = Positions(alteri,:)+round(ImgBox(IDsLDR(alteri),:)*1.1);
            Screen('DrawTexture', Win, itemTx{IDsLDR(alteri)}, [], Positions(alteri,:)+ImgBox(IDsLDR(alteri),:));
            Screen('FrameRect',Win,black,PresentBoxes(stacki,:),linewidth);
        end
        BoxL = PresentBoxes(1,:);
        BoxD = PresentBoxes(2,:);
        BoxR = PresentBoxes(3,:);
        
        t.options = Screen('Flip', Win, t.fix + ITI(trial) - ifi/2,1); % show stimuli
        
        KeyCode(KeyCode==1) = 0;
        while ~(KeyCode(bL) || KeyCode(bD) || KeyCode(bR)) && GetSecs < t.options + TPvalue(TimePressure) - ifi
            [~, t.choice, KeyCode, RTvariance] = KbCheck;
            Escape(KeyCode, escapeKey);
        end
        if ~(KeyCode(bL) || KeyCode(bD) || KeyCode(bR))
            message = 'Time Out!';
            DrawFormattedText(Win, double(message),'center','center',red,[],[],[],2);
            t.choice = NaN;
            RT = NaN;
            RTvariance = NaN;
            chosenPosition = NaN;
            chosenID = NaN;
            chosenItem = NaN;
            chosenValue = NaN;
            t.feedback = Screen('Flip', Win);
            WaitSecs(2.0);
            timeout = 1;
        else
            choice = [KeyCode(bL) KeyCode(bD) KeyCode(bR)];
            RT = t.choice - t.options;
            SlctBox = choice*[BoxL; BoxD; BoxR];
            [SlctCenterx,~] = RectCenter(SlctBox);
            chosenPosition = buttoncode{choice==1};
            chosenID = IDsLDR(SlctCenterx == Positions(:,1));
            chosenItem = find(IDs123 == chosenID);
            chosenValue = ValuesLDR(SlctCenterx == Positions(:,1));
            Screen('FrameRect',Win,red,SlctBox,linewidth);
            t.feedback = Screen('Flip', Win);
            %imgfile =  fullfile(scrn_dir,sprintf('ChoiceTrial_%03i.jpg',trial));
            %ImgArray = Screen('GetImage',Win);
            %imwrite(ImgArray,imgfile);
            WaitSecs(.5);
            timeout = 0;
        end
        fprintf(datafilepointer,'%i\t%s\t%i\t%s\t%2.1f\t%s\t%i\t%i\t%i\t%i\t%.1f\t%.1f\t%.1f\t%i\t%i\t%i\t%.1f\t%.1f\t%.1f\t%i\t%.3f\t%.3f\t%s\t%i\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n',...
            group, subID, trial, TPcode{TimePressure}, TPvalue(TimePressure), Vaguecode{DsnMtx.Vagueness(trial)}, DsnMtx.Vagueness(trial),...
        IDs123(1),IDs123(2),IDs123(3),Values123(1),Values123(2),Values123(3),...
        IDsLDR(1), IDsLDR(2), IDsLDR(3), ValuesLDR(1), ValuesLDR(2), ValuesLDR(3),...
        chosenItem, RT, RTvariance, chosenPosition, chosenID, chosenValue,...
            t.fix-t.main, t.options-t.main, t.choice-t.main, t.feedback-t.main, timeout);
        %% session break
        if any(trial == breakat-1)
            Screen('DrawTexture', Win, FixTx, [], FixationBox);
            t.inbreak = Screen('Flip', Win, t.feedback + .5 - ifi/2);
            WaitSecs(ITI(trial+1));
            message = sprintf('Take a short break...\nPress any key to continue when you are ready');
            DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
            Screen('Flip', Win);
            % Waiting for any key press
            t.enter = KbStrokeWait;
        end
    end
    copyfile(logtxt,DatPckg_dir);
    message = 'Task Complete!\nNow, relax\nWait for your gift!';
    DrawFormattedText(Win,double(message),'center','center', black,[],[],[],2);
    Screen('Flip', Win, t.feedback + .5);
    
    % loop until 'enter' key is pressed, controlled by experimenter
    KeyCode(KeyCode==1) = 0;
    while ~KeyCode(enter_button)
        [~, ~, KeyCode, ~] = KbCheck;
        Escape(KeyCode, escapeKey);
    end
    sca;
    ShowCursor;
    fclose('all');
    Priority(0);
catch
    % catch error: This is executed in case something goes wrong in the
    % 'try' part due to programming error etc.:
    % Do same cleanup as at the end of a regular session...
    sca; %Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end

end
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