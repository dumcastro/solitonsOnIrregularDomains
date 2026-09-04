function ax = wave_colorstrip(fig, pos, cmap, clim, gamma, label)
%WAVE_COLORSTRIP  Horizontal colour key drawn as its own axes.
%
%   AX = WAVE_COLORSTRIP(FIG,POS,CMAP,CLIM,GAMMA,LABEL) draws a horizontal
%   colour bar for CMAP at figure-normalised position POS = [x y w h].
%
%   Unlike COLORBAR this does not resize a neighbouring axes, so it can be
%   shared cleanly by a whole montage of panels.  When GAMMA is not 1 the
%   ticks sit at the stretched positions of round values of w and carry
%   those true values as labels, so the key stays honest about the
%   nonlinear colour mapping.
%
%   See also PLOT_SOLUTION, WAVE_COLORBAR, WAVE_COLORMAP.

if nargin < 5 || isempty(gamma), gamma = 1; end
if nargin < 6, label = ''; end
A = max(abs(clim));

n  = size(cmap, 1);
ax = axes('Parent', fig, 'Position', pos);
image(ax, linspace(-1, 1, n), [0 1], reshape(1:n, [1 n]));
colormap(ax, cmap);
set(ax, 'CLim', [1 n], 'YTick', [], 'TickDir', 'out', 'Box', 'on', ...
        'Layer', 'top', 'FontSize', 9);

frac = [-1 -0.5 -0.2 0 0.2 0.5 1];
vals = frac*A;
pp   = sign(frac).*abs(frac).^gamma;
[pp, iu] = unique(pp);
vals = vals(iu);
lab  = cell(1, numel(vals));
for i = 1:numel(vals), lab{i} = sprintf('%.3g', vals(i)); end
set(ax, 'XTick', pp, 'XTickLabel', lab);

if ~isempty(label)
    xlabel(ax, label, 'FontSize', 10);
end
