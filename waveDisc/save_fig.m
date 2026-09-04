function save_fig(fig, fname, dpi)
%SAVE_FIG  Write a figure to PNG, using the best method the release offers.
%
%   SAVE_FIG(FIG,FNAME) writes FIG to FNAME at 200 dpi.  EXPORTGRAPHICS is
%   used where it exists (R2020a+, tight cropping and true colour), with
%   PRINT as the fallback for older MATLAB and for Octave.
if nargin < 3 || isempty(dpi), dpi = 200; end
d = fileparts(fname);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
try
    exportgraphics(fig, fname, 'Resolution', dpi, 'BackgroundColor', 'white');
catch
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    print(fig, fname, '-dpng', sprintf('-r%d', dpi));
end
fprintf('    wrote %s\n', fname);
end
