function Color_main(SubID,SubName, op_JND)
% Set Default JND
if nargin > 2
    JND = op_JND;
else
    JND = 1;
end
%% import participant data
filename = sprintf('Results_%d_%s', SubID, SubName);
%% Design Matrix
rng('shuffle');
frame_rect = [0 0 16 9]; % the 16:9 outside frame
velo_control = 4.5; % velocity of control dot,
num_conds = 8; % number of conditions in constant stimuli method
velo_exp = linspace(velo_control - JND*7/4, velo_control + JND*7/4, num_conds);
%velo_exp = [2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25]; % velocity of test dot in deg/sec
trial_per_condition = 48; % repetition number for each condition
trial_num = trial_per_condition * length(velo_exp); % total trial number
trial_rest = linspace(0,trial_num,5); % time for rest
all_conds = genTrials(trial_per_condition/8,3,[length(velo_exp), 2, 2, 2]); % generate all conds
all_conds = all_conds(all_conds(:,1)==1,:);
all_conds = all_conds(randperm(size(all_conds,1)),:);
cond = all_conds(:,2); % generate the condition for each trial
correct_ans = velo_exp(cond) > velo_control; % the correct answer (< or >)
exp_cond = all_conds(:,3);
inf_cond = all_conds(:,4);
aud_cond = all_conds(:,5);
blue = [0,0,255]; red = [255,0,0]; green = [0,255,0];
colors = [{green}, {red}; {red}, {green}];
% in color differentiation: 1 - red is exp, green is control;
%                           2 - red is control, green is exp;
inf_color = blue; % for interference dot
% initialize data recording variables
res = zeros(1, trial_num); % direct response, corresponds to exp_cond
RTs = zeros(1, trial_num); % Reaction time for each trial

%% setups
% Keyboard Recording setup
KbName('UnifyKeyNames');
key_esc           =   KbName('ESCAPE');
key_space           =   KbName('space');
key_f = KbName('f'); key_j = KbName('j'); 
key_q = KbName('q'); key_p = KbName('p');
SubPressMapping = [key_f,key_j;key_j,key_f];
% odd_num: f-green; j-red; even_num: f-red; j-green;
key_red = SubPressMapping(mod(SubID,2)+1,1);
key_green = SubPressMapping(mod(SubID,2)+1,2);
RestrictKeysForKbCheck([key_esc,key_space,key_red,key_green,key_q,key_p]);

instruction1 = 'If green is faster, press J; if red is faster, press F.';
instruction2 = 'If green is faster, press F; if red is faster, press J.';
intro = eval(sprintf('instruction%.d;',mod(SubID,2)+1));
startwithspace = 'Press space to start.';
response_text = 'Response Now!';

try
    %% %%%%%%%%%%% Visual Setup %%%%%%%%%%%%%%%%%%%
    % Open up a screen
    Screen('Preference', 'SkipSyncTests', 1);
    ScreenNumber = max(Screen('Screens')); % count the screens
    AssertOpenGL;
    InitializeMatlabOpenGL;

    % Initiating visual display
    background_color = [0, 0, 0]; % black background
    [w, Rect] = Screen('OpenWindow',ScreenNumber,background_color);
    [cx, cy] = RectCenter(Rect); % get the center point
    white = WhiteIndex(w); black = BlackIndex(w);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [ResX, ResY] = Screen('WindowSize',w);
    [width, height]=Screen('DisplaySize', w);
    flipIntv = Screen('GetFlipInterval',w);
    slack = flipIntv / 2;
    RefreshRate = 60; % get frame rate (about 100 frame/sec)
    vdist = 50; % set visual distance as 50 cm
    pxlpdg = deg2pix(1,sqrt(width^2+height^2)/25.4,ResX,vdist,ResY/ResX); % calculate pixel per degree
    Screen(w,'TextSize', 30);
    Screen('TextFont',w,'Arial');
    Screen(w,'TextStyle',0); % 1bold text;0normal
    HideCursor;

    %%%%% locate the frame %%%%%
    scale = round((ResY/2)/9); % let the frame as half height as the screen
    scale_frame = ScaleRect(frame_rect, scale, scale); % scale the frame
    centered_frame = CenterRectOnPoint(scale_frame, cx, cy); % locate the frame at center
    frame_xlim = centered_frame([1,3]); % the horizontal boundary of frame
    frame_ylim = centered_frame([2,4]); % the vertical boundary of frame

    %%%%% initialize the dot stimuli %%%%%
    dot_rect = [0 0 20 20]; % with radius = 10  pixel
    %% %%%%%%%%%% Initiation audiotory display %%%%%%%%
    % Parameters for leading
    sr = 44100; % sampling rate
    InitializePsychSound(1);
    PsychPortAudio('Close');
    pahandle = PsychPortAudio('Open', [], [], [], sr, 2);
    PsychPortAudio('Volume',pahandle, 0.1);

    %%%%% initialize the auditory stimuli %%%%%
    dur = 10 * flipIntv; % length in sec (10 frames)
    freq = 1000; % frequency of the auditory stimulus
    gatedur = 2 * flipIntv; % length of gate
    tone = Soundgenerate(sr,freq,dur,gatedur);

    %% %%%%% calculate the trajectories and auditory stimulus %%%%%
    trajectories = cell(trial_num,1); % create a cell recording the trajectories
    sounds = cell(trial_num, 1); % create a cell recording the sounds

    for i = 1:trial_num
        %%% randomize the trajectories %%%
        n_change = 2; % each dot change its orientation for 3 times
        trial_dur = n_change * RefreshRate; % each trial last for 2 secs (120 frames)

        % each row indicates one dot
        change_occur = randi([30,60],2,n_change); % generate the time when orientation changes
        change_j = cumsum(change_occur, 2); % js at which orientation changes
        step_len = [change_occur, trial_dur - sum(change_occur, 2)]; % length for each step

        velo_ = [velo_control/RefreshRate; velo_exp(cond(i))/RefreshRate]; % velo in frame
        velo_ = velo_ .* pxlpdg; % for visible movement

        ori_change_lim = [60, 120]; % the range of orientation change in degree
        ori_ = zeros(2, n_change + 1);
        ori_(:,1) = randi(360,2,1);% the initial moving orientation

        % add two possible interference dots
        inf_change_occ = randi([30,60],2,n_change);
        inf_change_j = cumsum(inf_change_occ, 2);
        inf_step_len = [inf_change_occ, trial_dur - sum(inf_change_occ, 2)];
        inf_velo = [2.75 + rand() * 3.5;2.75 + rand() * 3.5];
        inf_velo_ = inf_velo / RefreshRate * pxlpdg;
        inf_ori_ = zeros(2, n_change + 1);
        inf_ori_(:,1) = randi(360,2,1);

        for j = 1:n_change
            %calculate the orientation for next step
            ori_(:,j+1) = ori_(:,j) + randi(ori_change_lim,2,1);
            inf_ori_(:,j+1) = inf_ori_(:,j) + randi(ori_change_lim,2,1);
        end

        dx = velo_ .* cos(ori_ .* pi ./ 180); % pixel change in x per frame
        dy = velo_ .* sin(ori_ .* pi ./ 180); % pixel change in y per frame

        inf_dx = inf_velo_ .* cos(inf_ori_ .* pi ./ 180);
        inf_dy = inf_velo_ .* sin(inf_ori_ .* pi ./ 180);

        %%% the initial position of dots %%%
        % prevent the out-of-boundary movement
        initial_xlim = [frame_xlim(1)+60,frame_xlim(2)-60];
        inf_xlim = [frame_xlim(1)+80,frame_xlim(2)-80];
        initial_ylim = [frame_ylim(1)+60,frame_ylim(2)-60];
        inf_ylim = [frame_ylim(1)+80,frame_ylim(2)-80];
        restrict_rect = [initial_xlim(1), initial_ylim(1), initial_xlim(2), initial_ylim(2)];

        %%% calculate the trajectory %%%
        %initialize
        traj_cont = zeros(trial_dur,2); traj_exp = zeros(trial_dur,2);
        traj_inf1 = zeros(trial_dur,2); traj_inf2 = zeros(trial_dur,2);
        traj_cont(1,:) = [randi(initial_xlim), randi(initial_ylim)];
        traj_exp(1,:) = [randi(initial_xlim), randi(initial_ylim)];
        traj_inf1(1,:) = [randi(inf_xlim), randi(inf_ylim)];
        traj_inf2(1,:) = [randi(inf_xlim), randi(inf_ylim)];
        k_cont = 0; k_exp = 0; k_inf1 = 0; k_inf2 = 0;% step for dots

        for j = 1:trial_dur-1
            % judge the changing point
            if ismember(j, change_j(1,:))||~k_cont
                k_cont = k_cont + 1;
                % check whether the next step would exceed the boundary
                expect_cont_x = traj_cont(j,1) + step_len(1,k_cont) * dx(1,k_cont);
                expect_cont_y = traj_cont(j,2) + step_len(1,k_cont) * dy(1,k_cont);
                while ~IsInRect(expect_cont_x, expect_cont_y, restrict_rect)
                    % if the expected trajectory exceeds the boundary
                    ori_(1,k_cont) = ori_(1,k_cont) - 15; %adjust the orientation
                    % calculate the expected trajectory again
                    dx = velo_ .* cos(ori_ .* pi ./ 180); % pixel change in x per frame
                    dy = velo_ .* sin(ori_ .* pi ./ 180); % pixel change in y per frame
                    expect_cont_x = traj_cont(j,1) + step_len(1,k_cont) * dx(1,k_cont);
                    expect_cont_y = traj_cont(j,2) + step_len(1,k_cont) * dy(1,k_cont);
                end
            end
            if ismember(j, change_j(2,:))||~k_exp
                k_exp = k_exp + 1;
                % check whether the next step would exceed the boundary
                expect_exp_x = traj_exp(j,1) + step_len(2,k_exp) * dx(2,k_exp);
                expect_exp_y = traj_exp(j,2) + step_len(2,k_exp) * dy(2,k_exp);
                while ~IsInRect(expect_exp_x, expect_exp_y, restrict_rect)
                    % if the expected trajectory exceeds the boundary
                    ori_(2,k_exp) = ori_(2,k_exp) - 15; %adjust the orientation
                    % calculate the expected trajectory again
                    dx = velo_ .* cos(ori_ .* pi ./ 180); % pixel change in x per frame
                    dy = velo_ .* sin(ori_ .* pi ./ 180); % pixel change in y per frame
                    expect_exp_x = traj_exp(j,1) + step_len(2,k_exp) * dx(2,k_exp);
                    expect_exp_y = traj_exp(j,2) + step_len(2,k_exp) * dy(2,k_exp);
                end
            end
            if ismember(j, inf_change_j(1,:))||~k_inf1
                k_inf1 = k_inf1 + 1;
                % check whether the next step would exceed the boundary
                expect_inf1_x = traj_inf1(j,1) + inf_step_len(1,k_inf1) * inf_dx(1,k_inf1);
                expect_inf1_y = traj_inf1(j,2) + inf_step_len(1,k_inf1) * inf_dy(1,k_inf1);
                while ~IsInRect(expect_inf1_x, expect_inf1_y, restrict_rect)
                    % if the expected trajectory exceeds the boundary
                    inf_ori_(1,k_inf1) = inf_ori_(1,k_inf1) - 5; %adjust the orientation
                    % calculate the expected trajectory again
                    inf_dx = inf_velo_ .* cos(inf_ori_ .* pi ./ 180); % pixel change in x per frame
                    inf_dy = inf_velo_ .* sin(inf_ori_ .* pi ./ 180); % pixel change in y per frame
                    expect_inf1_x = traj_inf1(j,1) + inf_step_len(1,k_inf1) * inf_dx(1,k_inf1);
                    expect_inf1_y = traj_inf1(j,2) + inf_step_len(1,k_inf1) * inf_dy(1,k_inf1);
                end
            end
            if ismember(j, inf_change_j(1,:))||~k_inf2
                k_inf2 = k_inf2 + 1;
                % check whether the next step would exceed the boundary
                expect_inf2_x = traj_inf2(j,1) + inf_step_len(2,k_inf2) * inf_dx(2,k_inf2);
                expect_inf2_y = traj_inf2(j,2) + inf_step_len(2,k_inf2) * inf_dy(2,k_inf2);
                while ~IsInRect(expect_inf2_x, expect_inf2_y, restrict_rect)
                    % if the expected trajectory exceeds the boundary
                    inf_ori_(2,k_inf2) = inf_ori_(2,k_inf2) - 15; %adjust the orientation
                    % calculate the expected trajectory again
                    inf_dx = inf_velo_ .* cos(inf_ori_ .* pi ./ 180); % pixel change in x per frame
                    inf_dy = inf_velo_ .* sin(inf_ori_ .* pi ./ 180); % pixel change in y per frame
                    expect_inf2_x = traj_inf2(j,1) + inf_step_len(2,k_inf2) * inf_dx(2,k_inf2);
                    expect_inf2_y = traj_inf2(j,2) + inf_step_len(2,k_inf2) * inf_dy(2,k_inf2);
                end
            end

            % iteration
            traj_cont(j+1,:) = traj_cont(j,:) + [dx(1,k_cont), dy(1, k_cont)];
            traj_exp(j+1,:) = traj_exp(j,:) + [dx(2,k_exp), dy(2,k_exp)];
            traj_inf1(j+1,:) = traj_inf1(j,:) + [inf_dx(1,k_inf1), inf_dy(1,k_inf1)];
            traj_inf2(j+1,:) = traj_inf2(j,:) + [inf_dx(2,k_inf2), inf_dy(2,k_inf2)];
        end

        %%% calculate the sound %%%
        sound_dur = round(trial_dur * flipIntv * sr); % duration for sound
        sound_ = zeros(1,sound_dur); % initialize the sound
        tone_length = length(tone); % the length of the tone (*sr)
        tone_start = round(change_j(1,:) .* flipIntv .* sr); % starting point
        tone_stop = tone_start + tone_length - 1; % stop point
        for j = 1:length(tone_start)
            sound_(tone_start(j):tone_stop(j)) = tone;
        end
        sounds{i} = sound_;

        %%% save for the display %%%
        trajectory.cont = traj_cont;
        trajectory.exp = traj_exp;
        trajectory.inf1 = traj_inf1;
        trajectory.inf2 = traj_inf2;
        trajectory.n_change = n_change;
        trajectory.change_j = change_j;
        trajectory.trial_dur = trial_dur;
        trajectory.velocity = velo_;
        trajectory.orientation = ori_;
        trajectories{i} = trajectory;
    end

    %% %%%%%%%%%% Pesentation %%%%%%%%%%%%
    drawTextAt(w,intro,ResX/2,ResY/4,white); drawTextAt(w,startwithspace,cx,cy,white); Screen('Flip',w);
    KbStrokeWait();
    exitFlag = false; %define exit flag

    %%% practice trials %%%
    prac_velo_inf = [2.75 + rand() * 3.5;2.75 + rand() * 3.5];
    prac_cond = randi(8,1,2); % the velo cond 
    prac_correct_ans = velo_exp(prac_cond) > velo_control;
    prac_inf_col = blue;
    prac_color_cond = (rand(1,2) > 0.5)+1; % the color cond
    prac1_cont_col = colors{prac_color_cond(1),1};prac1_exp_col = colors{prac_color_cond(1),2};
    prac2_cont_col = colors{prac_color_cond(2),1};prac2_exp_col = colors{prac_color_cond(2),2};
    prac1_exp_traj = Single_Traj(velo_exp(prac_cond(1)),RefreshRate,centered_frame,pxlpdg,n_change);
    prac1_con_traj = Single_Traj(velo_control,RefreshRate,centered_frame,pxlpdg,n_change);
    prac2_exp_traj = Single_Traj(velo_exp(prac_cond(2)),RefreshRate,centered_frame,pxlpdg,n_change);
    prac2_con_traj = Single_Traj(velo_control,RefreshRate,centered_frame,pxlpdg,n_change);
    prac2_inf1_traj = Single_Traj(prac_velo_inf(1),RefreshRate,centered_frame,pxlpdg,n_change);
    prac2_inf2_traj = Single_Traj(prac_velo_inf(2),RefreshRate,centered_frame,pxlpdg,n_change);
    prac1_sound = Single_Sound(RefreshRate,tone,flipIntv,sr,prac1_con_traj.change_j);
    prac2_sound = Single_Sound(RefreshRate,tone,flipIntv,sr,prac2_con_traj.change_j);
    in_prac1 = true;
    in_prac2 = true;
    while in_prac1
        drawTextAt(w,'Practice trial one (without interference)',cx,cy-100,white);
        t = Screen('Flip',w);
        Screen('Flip',w,t+1);
        PsychPortAudio('FillBuffer', pahandle, [prac1_sound;prac1_sound]);
        PsychPortAudio('Start', pahandle, 1,0,0);
        for j = 1:prac1_con_traj.trial_dur
            Screen('FrameRect', w, white, centered_frame);
            loc_cont = prac1_con_traj.traj(j,:); % location of control dot
            loc_exp = prac1_exp_traj.traj(j,:); % location of experimental dot
            cont_dot = CenterRectOnPoint(dot_rect, loc_cont(1), loc_cont(2));
            exp_dot = CenterRectOnPoint(dot_rect, loc_exp(1), loc_exp(2));
            %%% draw visual stimulus %%%
            Screen('FillOval', w, [prac1_cont_col', prac1_exp_col'], [cont_dot', exp_dot']);
            Screen('Flip',w);
        end
        drawTextAt(w,response_text,ResX/2,ResY/8,white);
        Screen('FrameRect', w, white, centered_frame);
        t_start = Screen('Flip',w);
        prac1_res = 0;
        while ~exitFlag && GetSecs()-t_start <= 3
            [keyIsDown, Secs, keyCode] = KbCheck;
            if keyCode(key_esc)
                exitFlag = true;
            end
            if keyCode(key_red)
                prac1_res = 1;
                prac1_RT = Secs - t_start; break;
            elseif keyCode(key_green)
                prac1_res = 2;
                prac1_RT = Secs - t_start; break;
            end
        end
        Screen('Flip',w);
        if exitFlag; PsychPortAudio('Close'); Screen('CloseAll'); break; end
        if ~prac1_res
            prac1_RT = -1;
            output = 'You did not respond in time!';
        else
            prac1_select_exp = (prac1_res == prac_color_cond(1));
            prac1_yesno = (prac1_select_exp == prac_correct_ans(1));
            if prac1_yesno
                yesno_output = 'You are right!';
            else
                yesno_output = 'You are wrong!';
            end
            output = sprintf([yesno_output,' Reaction Time: %.2d'],prac1_RT);
        end 
        drawTextAt(w,'End of Practice trial 1',cx,cy-100,white);
        drawTextAt(w,output,cx,cy,white);
        drawTextAt(w,'press p to try again, press q to proceed',cx,cy+100,white);
        Screen('Flip',w);
        while 1
            [keyIsDown,~,keyCode]=KbCheck;
            if keyIsDown && keyCode(key_q)
                in_prac1 = false;
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

    while in_prac2
        drawTextAt(w,'Practice trial two (with interference)',cx,cy-100,white);
        t = Screen('Flip',w);
        Screen('Flip',w,t+1);
        PsychPortAudio('FillBuffer', pahandle, [prac2_sound;prac2_sound]);
        PsychPortAudio('Start', pahandle, 1,0,0);
        for j = 1:prac2_con_traj.trial_dur
            Screen('FrameRect', w, white, centered_frame);
            loc_cont = prac2_con_traj.traj(j,:); % location of control dot
            loc_exp = prac2_exp_traj.traj(j,:); % location of experimental dot
            loc_inf1 = prac2_inf1_traj.traj(j,:);
            loc_inf2 = prac2_inf2_traj.traj(j,:);
            cont_dot = CenterRectOnPoint(dot_rect, loc_cont(1), loc_cont(2));
            exp_dot = CenterRectOnPoint(dot_rect, loc_exp(1), loc_exp(2));
            inf1_dot = CenterRectOnPoint(dot_rect, loc_inf1(1), loc_inf1(2));
            inf2_dot = CenterRectOnPoint(dot_rect, loc_inf2(1), loc_inf2(2));
            %%% draw visual stimulus %%%
            Screen('FillOval', w, [prac2_cont_col', prac2_exp_col', prac_inf_col', prac_inf_col'], [cont_dot', exp_dot', inf1_dot', inf2_dot']);
            Screen('Flip',w);
        end
        drawTextAt(w,response_text,ResX/2,ResY/8,white);
        Screen('FrameRect', w, white, centered_frame);
        t_start = Screen('Flip',w);
        prac2_res = 0;
        while ~exitFlag && GetSecs()-t_start <= 3
            [keyIsDown, Secs, keyCode] = KbCheck;
            if keyCode(key_esc)
                exitFlag = true;
            end
            if keyCode(key_red)
                prac2_res = 1;
                prac2_RT = Secs - t_start; break;
            elseif keyCode(key_green)
                prac2_res = 2;
                prac2_RT = Secs - t_start; break;
            end
        end
        Screen('Flip',w);
        if exitFlag; PsychPortAudio('Close'); Screen('CloseAll'); break; end
        if ~prac2_res
            prac2_RT = -1;
            output = 'You did not respond in time!';
        else
            prac2_select_exp = (prac2_res == prac_color_cond(2));
            prac2_yesno = (prac2_select_exp == prac_correct_ans(2));
            if prac2_yesno
                yesno_output = 'You are right!';
            else
                yesno_output = 'You are wrong!';
            end
            output = sprintf([yesno_output,' Reaction Time: %.2d'],prac2_RT);
        end 
        drawTextAt(w,'End of Practice trial 2',cx,cy-100,white);
        drawTextAt(w,output,cx,cy,white);
        drawTextAt(w,'press p to try again, press q to proceed',cx,cy+100,white);
        Screen('Flip',w);
        while 1
            [keyIsDown,~,keyCode]=KbCheck;
            if keyIsDown && keyCode(key_q)
                in_prac2 = false;
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

    drawTextAt(w,intro,cx,cy/2,white);
    drawTextAt(w,'Press space to start the experiment...',cx,cy,white);
    Screen('Flip',w);
    KbStrokeWait();
    Screen('Flip',w);
    
    for i = 1:trial_num
        %%% for color differentiation %%%
        whether_inf = inf_cond(i);
        whether_aud = aud_cond(i);
        cont_col = colors{exp_cond(i),1}; exp_col = colors{exp_cond(i),2};
        PsychPortAudio('FillBuffer', pahandle, [sounds{i};sounds{i}]);
        if whether_aud == 1
            PsychPortAudio('Start', pahandle, 1, 0, 1);
        end
        for j = 1:trajectories{i}.trial_dur
            Screen('FrameRect', w, white, centered_frame);
            loc_cont = trajectories{i}.cont(j,:); % location of control dot
            loc_exp = trajectories{i}.exp(j,:); % location of experimental dot
            loc_inf1 = trajectories{i}.inf1(j,:); % location of interference dot
            loc_inf2 = trajectories{i}.inf2(j,:);
            cont_dot = CenterRectOnPoint(dot_rect, loc_cont(1), loc_cont(2));
            exp_dot = CenterRectOnPoint(dot_rect, loc_exp(1), loc_exp(2));
            inf1_dot = CenterRectOnPoint(dot_rect, loc_inf1(1), loc_inf1(2));
            inf2_dot = CenterRectOnPoint(dot_rect, loc_inf2(1), loc_inf2(2));
            %%% draw visual stimulus %%%
            if whether_inf == 1
                Screen('FillOval', w, [cont_col', exp_col', inf_color',inf_color'], [cont_dot', exp_dot', inf1_dot', inf2_dot']);
            else
                Screen('FillOval', w, [cont_col', exp_col'], [cont_dot', exp_dot']);
            end
            Screen('Flip',w);
            
        end
        drawTextAt(w,response_text,ResX/2,ResY/8,white);
        Screen('FrameRect', w, white, centered_frame);
        t_start = Screen('Flip',w);
        while ~exitFlag && GetSecs()-t_start <= 3
            [keyIsDown, Secs, keyCode] = KbCheck;
            if keyCode(key_esc)
                exitFlag = true;
            end
            if keyCode(key_red)
                res(i) = 1;
                RTs(i) = Secs - t_start; break;
            elseif keyCode(key_green)
                res(i) = 2;
                RTs(i) = Secs - t_start; break;
            end
        end
        Screen('FrameRect', w, white, centered_frame);
        Screen('Flip',w);
        if exitFlag; PsychPortAudio('Close'); Screen('CloseAll'); break; end
        if ~res(i); RTs(i) = -1; end % record for non-responding
        select_exp = (res == exp_cond'); % whether participant selected the exp dot
        yesno = (select_exp == correct_ans); % whether participant is correct
        WaitSecs(1+rand/2); %ITI

        if ismember(i, trial_rest(2:end-1))
            drawTextAt(w,'Please take a rest.',cx,cy-100,white);
            drawTextAt(w,'If you are ready to continue, please press space',cx,cy+100,white);
            Screen('Flip',w);
            KbStrokeWait();
            Screen('Flip',w);
        end
    end
    
    Screen('Flip',w);
    WaitSecs(1);
    ListenChar(0);
    ShowCursor;
    PsychPortAudio('Close');
    Screen('CloseAll');
    save(['Color_data/',filename,'.mat'], 'SubID', 'SubName',...
        "all_conds","inf_cond","exp_cond","correct_ans", "trajectories", "sounds", "RTs", "res" ,"yesno");

catch error
    ShowCursor;
    ListenChar(0)
    PsychPortAudio('Close');
    Screen('CloseAll');
    save(['Color_data/error_' filename '.mat'], 'SubID', 'SubName',...
        "all_conds","inf_cond","exp_cond", "correct_ans", "trajectories", "sounds" , "RTs" , "res")
    rethrow(error);
end

