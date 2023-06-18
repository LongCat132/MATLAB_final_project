function smoothed_tone= Soundgenerate(sr,frequency,dur,gatedur)
% Input: 
%   sr:sampling frequency; frequency: of the tone;
%   dur: stimulus duration; gatedur: gating duration
% Output:
%   the sound sequence of gated tone
time=linspace(0,dur,sr*dur);
tone=sin(2*pi*frequency*time);
gate = cos(linspace(pi, 2*pi, sr*gatedur)); % cosine window
gate = (gate + 1) / 2;
offsetgate = fliplr(gate); 
sustain = ones(1, (length(tone)-2*length(gate)));
envelope = [gate, sustain, offsetgate]; % add gates to the beginning and the end
smoothed_tone = envelope .* tone;
end
