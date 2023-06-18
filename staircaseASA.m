  function [intnext, done] = staircaseASA(intlist, resplist, threshold, stepinit, stepstop)
%%
%
% Accelerated stochastic approximation staircase. 
% see Faes et al 2007 and Garcia-Perez 2011
%
% input:
%   intlist   = list of intensities, including current intensity
%   resplist  = list of responses, including response to current intensity (0 = miss, 1 = hit)
%   threshold = desired threshold
%   stepinit  = maximum initial step size (in units of intensity)
%
% optional argument:
%   stepstop = lower limit for step as stopping criterion, otherwise just
%   call it 30x.
%
% output:
%   intnext = next intensity
%   done    = converged (1) or not yet (0)


nn = length(intlist);      % number of intensities displayed so far (including current)
intcurr = intlist(end);    % current intensity being displayed
cc = stepinit;

% cc = stepinit / max(threshold, 1-threshold);

%% accelerated part, finding the number of response flips. 
nshift = 0;                     % number of shifts in response categories
respprev = resplist(1);
for ii = 2:nn
    respcurr = resplist(ii);
    if (respcurr ~= respprev)
        nshift = nshift + 1;
    end
    respprev = respcurr; 
end

%% calculate next intensity
respcurr = resplist(end);
if (nn <= 2)
    intnext = intcurr - (cc / nn) * (respcurr - threshold);
else
    intnext = intcurr - (cc / (nshift + 1)) * (respcurr - threshold);
end

if (nargin > 4)
    done = (abs(intnext - intcurr) < stepstop);
end
