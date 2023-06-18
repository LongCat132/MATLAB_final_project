%%%%%%% Pre_Block Testing Individual JND
function [JND, done] = Speed_Judge(SubID,SubName)
% default params
% SubID = 0;
% SubName = 'John';
%addpath("__self_func__")
rng('shuffle');
%% Design Matrix
frame_rect = [0 0 4 4]; % the 4:4 outside frame; two frames
velo_control = 4.5; % velocity of control dot,
velo_delta_lst = [2];
maxtrialnum = 50;
grey = [128, 128, 128];
% initialize data recording variables
Response = zeros(3, maxtrialnum); % direct response, corresponds to exp_cond
% c1: exp speed; c2: resp_key(1-left fast, 2-right fast); c3: accuracy
RTs = zeros(1, maxtrialnum); % Reaction time for each trial
accuracy = -1*zeros(1, maxtrialnum);
%% setups
% ASA set up
threshold = 0.75;
stepinit = 5;
stepstop = 0.1;
% Keyboard Recording setup
KbName('UnifyKeyNames');
key_esc          =   KbName('ESCAPE');
key_space          =   KbName('space');
key_f = KbName('f'); key_j = KbName('j');
key_q = KbName('q'); key_p = KbName('p');
% RestrictKeysForKbCheck([key_esc,key_space,key_f,key_j,key_q,key_p]);
instruction = 'If Left is faster, press F; if Right is faster, press J.';
startwithspace = 'Press space to start.';
response_text = 'Response Now!';

%%
try
    %% %%%%%%%%%%% Visual Setup %%%%%%%%%%%%%%%%%%%
    % Open up a screen
    Screen('Preference', 'SkipSyncTests', 1);
    ScreenNumber = max(Screen('Screens')); % count the screens
    AssertOpenGL;
    InitializeMatlabOpenGL;

    % Initiating visual display
    background_color = [0, 0, 0]; % black background
    [w, Rect] = Screen('OpenWindow',ScreenNumber,background_color, [],[],2);
    [cx, cy] = RectCenter(Rect); % get the center point
    white = WhiteIndex(w); black = BlackIndex(w);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [ResX, ResY] = Screen('WindowSize',w);
    [width, height]=Screen('DisplaySize', w);
    flipIntv = Screen('GetFlipInterval',w);
    slack = flipIntv / 2;
    RefreshRate = 60; % get frame rate (about 60 frame/sec)
    vdist = 50; % set visual distance as 50 cm
    pxlpdg = deg2pix(1,sqrt(width^2+height^2)/25.4,ResX,vdist,ResY/ResX); % calculate pixel per degree
    Screen(w,'TextSize', 30);
    Screen('TextFont',w,'Arial');
    Screen(w,'TextStyle',0); % 1bold text;0normal
    cross = MakeCross(w,0.05,0.5,pxlpdg,[128 128 128],[0 0 0]);
    blue_cross = MakeCross(w,0.05,0.5,pxlpdg,[0 0 255],[0 0 0]);
    HideCursor;

    %%%%% locate the frame %%%%%
    scale = round(ResY/8); % scale the frame
    scale_frame = ScaleRect(frame_rect, scale, scale); % scale the frame
    Left_frame = CenterRectOnPoint(scale_frame, round(cx*0.6), cy); % locate the frame at center
    Right_frame = CenterRectOnPoint(scale_frame, round(cx*1.4), cy);
    Frames = {Left_frame, Right_frame; Right_frame, Left_frame};
    % 1- Left exp; 2- Right exp;
    %%%%% initialize the dot stimuli %%%%%
    dot_rect = [0 0 20 20]; % with radius = 10  pixel

    %% Presentation
    drawTextAt(w,['Welcome to our experiment, ',SubName,'!'],ResX/2, ResY/4, white);
    drawTextAt(w,instruction,ResX/2, ResY/2, white);
    drawTextAt(w,startwithspace,ResX/2, 5*ResY/8, white);
    Screen('Flip',w);
    keyIsDown = 0;
    while 1
        [keyIsDown,~,keyCode]=KbCheck;
        if keyIsDown && keyCode(key_space)
            break;
        elseif keyIsDown && keyCode(key_esc)
            Screen('CloseAll');
            ShowCursor;
            break;
        end
    end
    Screen('FrameRect', w, ones(3,2)*255, [Left_frame',Right_frame'],ones(1,2)*5);
    Screen('DrawTexture', w, cross);
    Screen('Flip',w);
    WaitSecs(1);


    %%%%%%% Presentation %%%%%%%%%
    % continue or again?
    prac_again = 1;

    while prac_again
        %% Loop through trial
        delta_change = [-1,1];
        for iTrial = [1:maxtrialnum]
            % presentation
            conorexp = randi(2);
            expfast = randi(2)-1;
            velo_exp = abs(delta_change(expfast+1)*velo_delta_lst(iTrial) + velo_control);
            if velo_exp > 6.5
                velo_exp = 6;
            end
            exp_frame = Frames{conorexp,1}; cont_frame = Frames{conorexp,2};
            traj_cont = Single_Traj(velo_control,RefreshRate,cont_frame,pxlpdg);
            traj_exp = Single_Traj(velo_exp,RefreshRate,exp_frame,pxlpdg);
            for j = 1:179
                Screen('FrameRect', w, ones(3,2)*255, [Left_frame',Right_frame'],ones(1,2)*5);
                loc_cont = traj_cont.traj(j,:); % location of control dot
                loc_exp = traj_exp.traj(j,:); % location of experimental dot
                cont_dot = CenterRectOnPoint(dot_rect, loc_cont(1), loc_cont(2));
                exp_dot = CenterRectOnPoint(dot_rect, loc_exp(1), loc_exp(2));
                %%% draw visual stimulus %%%
                Screen('FillOval', w, [grey', grey'], [cont_dot', exp_dot']);
                Screen('DrawTexture', w, cross);
                Screen('Flip',w);
            end

            % Collect Response
            StartTime = GetSecs;
            Screen('DrawTexture', w, blue_cross);
            Screen('FrameRect', w, ones(3,2)*255, [Left_frame',Right_frame'],ones(1,2)*5);
            Screen('Flip',w);

            Response(iTrial,1) = velo_exp;
            Response(iTrial,2) = conorexp;

            keyIsDown = 0;
            while 1
                [keyIsDown,~,keyCode]=KbCheck;
                if keyIsDown && keyCode(key_f)
                    RTs(iTrial) = GetSecs - StartTime;
                    Response(iTrial,3) = 1;
                    break;
                elseif keyIsDown && keyCode(key_j)
                    RTs(iTrial) = GetSecs - StartTime;
                    Response(iTrial,3) = 2;
                    break;
                elseif keyIsDown && keyCode(key_esc)
                    Screen('CloseAll');
                    ShowCursor;
                    break;
                end
            end

            Screen('DrawTexture', w, cross);
            Screen('FrameRect', w, ones(3,2)*255, [Left_frame',Right_frame'],ones(1,2)*5);
            Screen('Flip',w);

            % calculate accuracy
            choose_exp = (Response(iTrial,3) == conorexp);
            accuracy(iTrial) = (expfast&&choose_exp)||(~expfast&&~choose_exp);

            % Update the next exp speed
            len = length(velo_delta_lst);
            Speed_ASA = velo_delta_lst(max(1,len-29):len);
            ResponseSeq_ASA = accuracy(max(1,len-29):len);
            intlst = 10*log10(Speed_ASA);
            [delta_SpeedNext,done] =...
                staircaseASA(intlst, ResponseSeq_ASA, threshold, stepinit, stepstop);
            delta_SpeedNext = 10.^(delta_SpeedNext/10);
            if delta_SpeedNext > velo_control
                delta_SpeedNext = velo_control - 0.02;
            end
            velo_delta_lst = [velo_delta_lst,delta_SpeedNext];
            % ITI
            WaitSecs(1+0.5*rand);
        end

        % calculate individual discrimination level
        JND = mean(Speed_ASA);

        drawTextAt(w,'Congrats! You have finished the PreBlock!', ResX/2, ResY/4, white);
        drawTextAt(w,['Your JND is', num2str(round(JND,3)),', continue (q) or prac again (p)?'], ResX/2, ResY/2, white);
        Screen('Flip',w);

        keyIsDown = 0;
        while 1
            [keyIsDown,~,keyCode]=KbCheck;
            if keyIsDown && keyCode(key_q)
                prac_again = 0;
                break;
            elseif keyIsDown && keyCode(key_p)
                break;
            elseif keyIsDown && keyCode(key_esc)
                Screen('CloseAll');
                ShowCursor;
                break;
            end
        end
    end

    %% Save the data
    filename = ['pre_limit_',num2str(SubID),SubName];
    save(['Color_data/pre_limits/',filename,'.mat'], 'SubID', 'SubName',...
       "Response","RTs","velo_delta_lst","JND","done");
    Screen('CloseAll');
    ListenChar(0);
catch error
    Screen('CloseAll');
    filename = ['pre_limit_error_',num2str(SubID),SubName];
    save(['Color_data/pre_limits/',filename,'.mat'], 'SubID', 'SubName',...
       "Response","RTs","velo_delta_lst");
    ListenChar(0);
    ShowCursor;
    rethrow(error);
end

