


function sound = Single_Sound(RefreshRate,tone,flipIntv,sr,change_j)
trial_dur = length(change_j) * RefreshRate; % each trial last for 3 secs (180 frames)
sound_dur = round(trial_dur * flipIntv * sr); % duration for sound
sound = zeros(1,sound_dur); % initialize the sound
tone_length = length(tone); % the length of the tone (*sr)
tone_start = round(change_j .* flipIntv .* sr); % starting point
tone_stop = tone_start + tone_length - 1; % stop point
for j = 1:length(tone_start)
    sound(tone_start(j):tone_stop(j)) = tone;
end














