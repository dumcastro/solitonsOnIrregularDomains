function sol = demo_wave_polygon(varargin)
%DEMO_WAVE_POLYGON  Phase 2 demonstration: a wave on a conformally mapped domain.
%
%   DEMO_WAVE_POLYGON() solves the wave equation on the image of the unit
%   disc under the analytic conformal map z = zeta + 0.22 zeta^4, a rounded
%   triangle (the perturbation e^{3 i theta} on the unit circle gives the
%   three lobes), and writes 'fig8_polygon.png' to the output folder.  The
%   map is analytic, so the demonstration needs no Schwarz-Christoffel
%   toolbox; pass 'Vertices' to run the same pipeline on a true polygon
%   once the toolbox is on the path.
%
%   The pulse here is deliberately wider than the one on the disc.  Phase 2
%   assembles a dense weighted mass matrix, at O(K^2) quadrature work, so
%   the affordable basis is a few hundred modes rather than a few thousand,
%   and the initial data has to be resolvable within it.
%
%   Options: 'OutDir' ('figures'), 'T' (5), 'LambdaMax' (26), 'sigma' (0.3),
%            'MapFun'/'MapDiff', 'Vertices'.
%
%   The figure shows the conformal grid, the first few eigenfunctions of
%   the mapped domain, and snapshots of the evolution, all drawn in the
%   physical plane.
%
%   See also WAVE_SOLVER_POLYGON, TEST_WAVE_SOLVER_POLYGON, DEMO_WAVE_DISC.

p = struct('OutDir', 'figures', 'T', 5, 'LambdaMax', 26, 'sigma', 0.30, ...
           'MapFun', [], 'MapDiff', [], 'Vertices', []);
fn = fieldnames(p);
for i = 1:2:numel(varargin)
    ix = find(strcmpi(fn, varargin{i}), 1);
    if isempty(ix)
        error('demo_wave_polygon:unknownOption', 'Unknown option ''%s''.', varargin{i});
    end
    p.(fn{ix}) = varargin{i+1};
end
if isempty(p.MapFun) && isempty(p.Vertices)
    p.MapFun  = @(z) z + 0.22*z.^4;
    p.MapDiff = @(z) 1 + 0.88*z.^3;
end
if ~exist(p.OutDir, 'dir'), mkdir(p.OutDir); end

fprintf('\n================ demo_wave_polygon ================\n');
if isempty(p.Vertices)
    sol = wave_solver_polygon('MapFun', p.MapFun, 'MapDiff', p.MapDiff, ...
                              'T', p.T, 'nT', 161, 'LambdaMax', p.LambdaMax, ...
                              'x0', 0.34, 'y0', 0.10, 'sigma', p.sigma);
else
    sol = wave_solver_polygon('Vertices', p.Vertices, 'T', p.T, 'nT', 161, ...
                              'LambdaMax', p.LambdaMax, 'sigma', p.sigma);
end

f  = sol.map;
cl = sol.Wdisp*[-1 1];
g  = 0.6;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 1450 800]);
set(fig, 'PaperPositionMode', 'auto');

% --- (a) the conformal grid --------------------------------------------
ax = axes('Parent', fig, 'Position', [0.035 0.56 0.20 0.34]);
hold(ax, 'on');
for rr = 0.1:0.1:1
    zz = f(rr*exp(1i*linspace(0, 2*pi, 400)));
    plot(ax, real(zz), imag(zz), 'Color', [0.55 0.70 0.85], 'LineWidth', 0.7);
end
for tt = 0:pi/12:2*pi-1e-9
    zz = f(linspace(0, 1, 120)*exp(1i*tt));
    plot(ax, real(zz), imag(zz), 'Color', [0.93 0.72 0.45], 'LineWidth', 0.7);
end
plot(ax, sol.boundary(:, 1), sol.boundary(:, 2), 'Color', [0.12 0.12 0.16], ...
     'LineWidth', 1.4);
hold(ax, 'off'); axis(ax, 'equal'); axis(ax, 'off');
title(ax, '(a) the conformal grid', 'FontWeight', 'normal');

% --- (b) the first three eigenfunctions of the mapped domain -----------
Pp = mode_fields(sol, 1:3);
for k = 1:3
    ax = axes('Parent', fig, 'Position', [0.035 + k*0.205, 0.56, 0.20, 0.34]);
    Wk = Pp(:, :, k); a = max(abs(Wk(:)));
    hs = surf(ax, sol.X, sol.Y, zeros(size(Wk)), Wk);
    set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(ax, 2); axis(ax, 'equal'); axis(ax, 'off');
    colormap(ax, wave_colormap('wave')); set(ax, 'CLim', a*[-1 1]);
    hold(ax, 'on');
    plot3(ax, sol.boundary(:, 1), sol.boundary(:, 2), 1e-3*ones(size(sol.boundary, 1), 1), ...
          'Color', [0.12 0.12 0.16], 'LineWidth', 1.1);
    hold(ax, 'off');
    title(ax, sprintf('(b%d) eigenfunction \\mu_%d = %.4f', k, k, sol.mu(k)), ...
          'FontWeight', 'normal');
end

% --- (c) the evolution --------------------------------------------------
tv = p.T*[0 0.12 0.26 0.45];
for k = 1:4
    ax = axes('Parent', fig, 'Position', [0.035 + (k-1)*0.205, 0.115, 0.20, 0.34]);
    plot_solution(sol, tv(k), 'Axes', ax, 'CLim', cl, 'Gamma', g, ...
                  'Title', sprintf('(c%d) t = %.2f', k, tv(k)));
end
wave_colorstrip(fig, [0.36 0.045 0.28 0.020], wave_colormap('wave'), cl, g, ...
                'displacement  w');

ax = axes('Parent', fig, 'Position', [0 0.955 1 0.04], 'Visible', 'off');
text('Parent', ax, 'Units', 'normalized', 'Position', [0.5 0.3], ...
     'String', ['Phase 2: the wave equation on a conformally mapped domain' ...
                '   (weighted eigenproblem -\Delta\phi = \mu^2|f''|^2\phi)'], ...
     'HorizontalAlignment', 'center', 'FontSize', 13);

save_fig(fig, fullfile(p.OutDir, 'fig8_polygon.png'));
close(fig);
fprintf('================ done ================\n\n');
end

% ======================================================================
function P = mode_fields(sol, ks)
%MODE_FIELDS  The first eigenfunctions of the mapped domain, on the plot grid.
E   = sol.E;
rp  = sol.rp;
thc = [sol.thp, 2*pi];
P   = zeros(numel(rp), numel(thc), numel(ks));
for i = 1:numel(ks)
    v = sol.V(:, ks(i));
    W = zeros(numel(rp), numel(thc));
    for j = 1:E.K
        if v(j) == 0, continue, end
        B = besselj(E.m(j), E.lambda(j)*rp);
        if E.par(j) > 0, A = B*cos(E.m(j)*thc); else, A = B*sin(E.m(j)*thc); end
        W = W + v(j)*A;
    end
    P(:, :, i) = W;
end
end
