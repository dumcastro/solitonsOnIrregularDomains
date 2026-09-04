function cmap = wave_colormap(name, n)
%WAVE_COLORMAP  Vivid colormaps for wave fields.
%
%   CMAP = WAVE_COLORMAP            256-colour diverging 'wave' map.
%   CMAP = WAVE_COLORMAP(NAME)      NAME is one of
%                                     'wave'     diverging navy-white-crimson
%                                     'aurora'   diverging teal-white-magenta
%                                     'ember'    sequential black-red-gold
%                                     'abyss'    sequential navy-cyan-white
%   CMAP = WAVE_COLORMAP(NAME,N)    N colours (default 256).
%
%   The diverging maps are centred on their middle entry, so pair them with
%   a symmetric colour axis, CAXIS([-A A]), to keep zero neutral.

if nargin < 1 || isempty(name), name = 'wave'; end
if nargin < 2 || isempty(n),    n    = 256;    end

switch lower(name)
    case 'wave'
        anchors = [ 0.031 0.098 0.290
                    0.086 0.271 0.580
                    0.157 0.502 0.796
                    0.388 0.729 0.878
                    0.729 0.898 0.929
                    0.984 0.984 0.973
                    0.996 0.878 0.545
                    0.973 0.667 0.271
                    0.898 0.400 0.169
                    0.741 0.129 0.129
                    0.451 0.020 0.086 ];
    case 'aurora'
        anchors = [ 0.000 0.180 0.204
                    0.000 0.396 0.435
                    0.098 0.631 0.616
                    0.502 0.831 0.749
                    0.851 0.945 0.906
                    1.000 1.000 1.000
                    0.976 0.855 0.925
                    0.918 0.600 0.816
                    0.780 0.310 0.667
                    0.514 0.129 0.478
                    0.243 0.043 0.267 ];
    case 'ember'
        anchors = [ 0.020 0.020 0.078
                    0.161 0.047 0.243
                    0.400 0.063 0.322
                    0.643 0.129 0.290
                    0.843 0.278 0.196
                    0.949 0.478 0.129
                    0.988 0.686 0.196
                    0.996 0.867 0.451
                    1.000 0.976 0.808 ];
    case 'abyss'
        anchors = [ 0.012 0.031 0.157
                    0.043 0.145 0.365
                    0.055 0.294 0.545
                    0.098 0.463 0.667
                    0.243 0.631 0.741
                    0.482 0.780 0.800
                    0.729 0.890 0.867
                    0.929 0.965 0.949 ];
    otherwise
        error('wave_colormap:badName', 'Unknown colormap ''%s''.', name);
end

na = size(anchors, 1);
t  = linspace(0, 1, na)';
tq = linspace(0, 1, n)';
cmap = zeros(n, 3);
for k = 1:3
    cmap(:, k) = interp1(t, anchors(:, k), tq, 'pchip');
end
cmap = min(max(cmap, 0), 1);
