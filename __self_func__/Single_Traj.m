function trajectory = Single_Traj(velo, RefreshRate, frame, pxlpdg, changes)
if nargin > 4
    n_change = changes;
    trial_secs = n_change;
else
    n_change = randi(3)+1;
    trial_secs = 3;
end

% velo = 5;
% RefreshRate = 60;
% frame = [200,200,800,800];
% pxlpdg = 42;

frame_xlim = frame([1,3]); % the horizontal boundary of frame
frame_ylim = frame([2,4]); % the vertical boundary of frame

%%% randomize the trajectories %%%
trial_dur = trial_secs * RefreshRate; % each trial last for 3 secs (180 frames)

% each row indicates one dot
change_occur = randi([30,60],1,n_change); % generate the time when orientation changes
change_j = cumsum(change_occur, 2); % js at which orientation changes
step_len = [change_occur, trial_dur - sum(change_occur, 2)]; % length for each step

velo_ = velo / RefreshRate * pxlpdg; % velo in frame

ori_change_lim = [60, 120]; % the range of orientation change in degree
ori_ = zeros(1, n_change + 1);
ori_(:,1) = randi(360);% the initial moving orientation

for j = 1:n_change
    %calculate the orientation for next step
    ori_(:,j+1) = ori_(:,j) + randi(ori_change_lim);
end

dx = velo_ .* cos(ori_ .* pi ./ 180); % pixel change in x per frame
dy = velo_ .* sin(ori_ .* pi ./ 180); % pixel change in y per frame

%%% the initial position of dots %%%
% prevent the out-of-boundary movement
initial_xlim = [frame_xlim(1)+60,frame_xlim(2)-60];
initial_ylim = [frame_ylim(1)+60,frame_ylim(2)-60];
restrict_rect = [initial_xlim(1), initial_ylim(1), initial_xlim(2), initial_ylim(2)];

%%% calculate the trajectory %%%
%initialize
traj = zeros(trial_dur,2);
traj(1,:) = [randi(initial_xlim,1), randi(initial_ylim,1)];
k = 0; % step for dots

for j = 1:trial_dur-1
    % judge the changing point
    if ismember(j, change_j)||~k
        k = k + 1;
        % check whether the next step would exceed the boundary
        expect_x = traj(j,1) + step_len(k) * dx(k);
        expect_y = traj(j,2) + step_len(k) * dy(k);
        while ~IsInRect(expect_x, expect_y, restrict_rect)
            % if the expected trajectory exceeds the boundary
            ori_(k) = ori_(k) - 15; %adjust the orientation
            % calculate the expected trajectory again
            dx = velo_ .* cos(ori_ .* pi ./ 180); % pixel change in x per frame
            dy = velo_ .* sin(ori_ .* pi ./ 180); % pixel change in y per frame
            expect_x = traj(j,1) + step_len(k) * dx(k);
            expect_y = traj(j,2) + step_len(k) * dy(k);
        end
    end
    % iteration
    traj(j+1,:) = traj(j,:) + [dx(1,k), dy(1, k)];
    trajectory.traj = traj;
    trajectory.n_change = n_change;
    trajectory.change_j = change_j;
    trajectory.trial_dur = trial_dur;
    trajectory.velocity = velo_;
    trajectory.orientation = ori_;
end

