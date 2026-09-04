function cb = wave_colorbar(ax, clim, gamma, label)
%WAVE_COLORBAR  Colourbar labelled in true field units.
%
%   CB = WAVE_COLORBAR(AX,CLIM,GAMMA) attaches a colourbar to AX for a plot
%   drawn by PLOT_SOLUTION with colour limits CLIM and colour stretch GAMMA.
%   When GAMMA is 1 this is an ordinary colourbar.  When GAMMA < 1 the
%   colours encode sign(w)|w/A|^gamma, so the ticks are placed at the
%   stretched positions of a set of round values of w and labelled with
%   those true values - the bar stays honest about what the colours mean.
%
%   WAVE_COLORBAR(AX,CLIM,GAMMA,LABEL) also sets the colourbar label.
%
%   See also PLOT_SOLUTION, WAVE_COLORMAP, SAFE_COLORBAR.

if nargin < 3 || isempty(gamma), gamma = 1;  end
if nargin < 4, label = ''; end
A = max(abs(clim));

cb = safe_colorbar(ax);

frac = [-1 -0.5 -0.2 -0.05 0 0.05 0.2 0.5 1];
vals = frac*A;
pos  = sign(vals).*(abs(vals)/A).^gamma*A;
[pos, iu] = unique(pos);
vals = vals(iu);

lab = cell(1, numel(vals));
for i = 1:numel(vals)
    lab{i} = sprintf('%.3g', vals(i));
end
try
    set(cb, 'Ticks', pos, 'TickLabels', lab);
catch
    set(cb, 'YTick', pos, 'YTickLabel', lab);   % older releases / Octave
end
if ~isempty(label)
    try, set(get(cb, 'Label'), 'String', label); catch, ylabel(cb, label); end
end
