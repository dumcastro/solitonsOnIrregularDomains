function h = plot_solution(sol, t, varargin)
%PLOT_SOLUTION  Render the disc wave field at a given time.
%
%   PLOT_SOLUTION(SOL,T) draws w(.,.,T) as a filled colour map on the domain.
%   SOL is the struct returned by WAVE_SOLVER_DISC or, equally, the one
%   returned by WAVE_SOLVER_POLYGON - the outline and the axis limits are
%   taken from SOL.boundary when it is present and from the disc radius
%   otherwise.  If T coincides with one of SOL.t the stored field is
%   reused; otherwise it is synthesised on the fly, so any T in the
%   interval may be asked for.
%
%   PLOT_SOLUTION(SOL,T,'Name',Value,...) accepts:
%     'Style'      'flat' (default) | 'surface' | 'contour'
%     'Axes'       axes handle to draw into                       (gca)
%     'CLim'       colour limits.  The default is symmetric at SOL.Wdisp,
%                  a robust amplitude for the whole run, so every frame of
%                  an animation shares one scale and zero stays neutral.
%     'Gamma'      signed power-law colour stretch, sign(w)|w/A|^gamma with
%                  A = max(CLim).  GAMMA < 1 lifts the weak ripples into
%                  view without changing the sign or the ordering of the
%                  data; the colourbar is then labelled with the true
%                  values of w at their stretched positions.          (1)
%     'Colormap'   name for WAVE_COLORMAP, or an N-by-3 matrix    ('wave')
%     'ZScale'     vertical exaggeration for 'surface'                (1)
%     'ZLim'       z limits for 'surface'.  Default autoscales to this
%                  frame; pass a fixed pair to keep an animation steady.
%     'Boundary'   draw the rim r = R                              (true)
%     'Lighting'   add a headlight to 'surface'.  Off by default, as it
%                  desaturates the colormap.                      (false)
%     'Title'      title string; '' suppresses it                  (auto)
%     'NContour'   number of levels for 'contour'                    (24)
%
%   H is a struct of the graphics handles that were created; H.clim and
%   H.gamma record the mapping actually used, which WAVE_COLORBAR needs.
%
%   Example
%     sol = wave_solver_disc('T', 6);
%     figure; plot_solution(sol, 2.4, 'Style', 'surface', 'ZScale', 1.5);
%
%   See also WAVE_SOLVER_DISC, WAVE_COLORMAP, WAVE_COLORBAR, DEMO_WAVE_DISC.

p = struct('Style', 'flat', 'Axes', [], 'CLim', [], 'Gamma', 1, ...
           'Colormap', 'wave', 'ZScale', 1, 'ZLim', [], 'Boundary', true, ...
           'Lighting', false, 'Title', [], 'NContour', 24);
if mod(numel(varargin), 2) ~= 0
    error('plot_solution:badArgs', 'Options must be Name/Value pairs.');
end
f = fieldnames(p);
for i = 1:2:numel(varargin)
    ix = find(strcmpi(f, varargin{i}), 1);
    if isempty(ix)
        error('plot_solution:unknownOption', 'Unknown option ''%s''.', varargin{i});
    end
    p.(f{ix}) = varargin{i+1};
end

W = field_at(sol, t);

if isempty(p.CLim)
    A = sol.Wdisp;
    if ~(A > 0), A = 1; end
    p.CLim = [-A A];
end
A = max(abs(p.CLim));

% Colour variable: the field itself, or its signed power-law stretch.
if p.Gamma == 1
    C = W;
else
    C = sign(W).*(min(abs(W)/A, 1).^p.Gamma)*A;
end

if ischar(p.Colormap), cmap = wave_colormap(p.Colormap); else, cmap = p.Colormap; end
if isempty(p.Axes), p.Axes = gca; end
ax = p.Axes;

X = sol.X; Y = sol.Y;
% The disc solver stores a radius; the polygon solver of Phase 2 stores the
% image of the unit circle instead.  Either way the outline and the axis
% limits come from the same two variables below.
if isfield(sol, 'boundary') && ~isempty(sol.boundary)
    bx = sol.boundary(:, 1); by = sol.boundary(:, 2);
else
    tb = linspace(0, 2*pi, 512)';
    bx = sol.opt.R*cos(tb); by = sol.opt.R*sin(tb);
end
pad = 0.04*max(max(bx) - min(bx), max(by) - min(by));
xl  = [min(bx) - pad, max(bx) + pad];
yl  = [min(by) - pad, max(by) + pad];
h = struct('clim', p.CLim, 'gamma', p.Gamma);

switch lower(p.Style)

    case 'flat'
        h.surf = surf(ax, X, Y, zeros(size(W)), C);
        set(h.surf, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(ax, 2);
        axis(ax, 'equal');
        set(ax, 'XLim', xl, 'YLim', yl);
        axis(ax, 'off');

    case 'surface'
        Z = p.ZScale*W;
        h.surf = surf(ax, X, Y, Z, C);
        set(h.surf, 'EdgeColor', 'none', 'FaceColor', 'interp');
        axis(ax, 'equal');
        if isempty(p.ZLim)
            zm = 1.15*max(max(abs(Z(:))), eps);
            p.ZLim = [-zm zm];
        end
        set(ax, 'XLim', xl, 'YLim', yl, 'ZLim', p.ZLim);
        view(ax, -37.5, 38);
        grid(ax, 'on');
        set(ax, 'Box', 'off');
        if p.Lighting
            try %#ok<TRYNC>
                camlight(ax, 'headlight'); lighting(ax, 'gouraud');
                material(h.surf, [0.55 0.55 0.12 8 0.4]);
            end
        end

    case 'contour'
        lv = linspace(p.CLim(1), p.CLim(2), p.NContour);
        [~, h.contour] = contourf(ax, X, Y, C, lv);
        set(h.contour, 'LineStyle', 'none');
        axis(ax, 'equal');
        set(ax, 'XLim', xl, 'YLim', yl);
        axis(ax, 'off');

    otherwise
        error('plot_solution:badStyle', ...
            'Style must be ''flat'', ''surface'' or ''contour''.');
end

colormap(ax, cmap);
set(ax, 'CLim', p.CLim);

if p.Boundary
    hold(ax, 'on');
    if strcmpi(p.Style, 'surface')
        h.rim = plot3(ax, bx, by, zeros(size(bx)));
    else
        h.rim = plot3(ax, bx, by, 1e-3*ones(size(bx)));
    end
    set(h.rim, 'Color', [0.12 0.12 0.16], 'LineWidth', 1.2);
    hold(ax, 'off');
end

if isempty(p.Title), p.Title = sprintf('t = %.2f', t); end
if ~isempty(p.Title)
    h.title = title(ax, p.Title);
    set(h.title, 'FontWeight', 'normal');
end
end

% ======================================================================
function W = field_at(sol, t)
%FIELD_AT  Stored snapshot when T is an output time, synthesised otherwise.
if ~isempty(sol.W)
    [dt, i] = min(abs(sol.t - t));
    if dt < 1e-12*max(1, abs(t))
        W = sol.W(:, :, i);
        return
    end
end
W = sol.field(t);
end
