function cb = safe_colorbar(ax)
%SAFE_COLORBAR  Attach a colourbar to AX across MATLAB and Octave releases.
%
%   COLORBAR(AX) is the modern spelling, COLORBAR('peer',AX) the old one,
%   and which of the two an installation accepts has changed more than once.
%   Try them in order and fall back to a bare COLORBAR on the current axes.
try
    cb = colorbar(ax);
catch
    try
        cb = colorbar('peer', ax);
    catch
        axes(ax); %#ok<LAXES>
        cb = colorbar;
    end
end
