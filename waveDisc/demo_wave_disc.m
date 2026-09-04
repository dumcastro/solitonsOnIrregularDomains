function sol = demo_wave_disc(varargin)
%DEMO_WAVE_DISC  Full visual demonstration of the disc wave solver.
%
%   DEMO_WAVE_DISC() solves the wave equation on the unit disc for a
%   Gaussian pulse released off centre and writes a set of figures and an
%   animation to the folder 'figures':
%
%     fig1_evolution.png     eight snapshots of the propagating front
%     fig2_surfaces.png      the same field as 3D relief
%     fig3_eigenmodes.png    gallery of the lowest Dirichlet eigenmodes
%     fig4_spectrum.png      where the initial data sits in the spectrum
%     fig5_verification.png  convergence, conservation and accuracy
%     fig6_spacetime.png     space-time diagrams and the characteristics
%     fig7_dirichlet_vs_neumann.png   the boundary-condition toggle
%     wave_disc.gif          animation of the evolution
%
%   DEMO_WAVE_DISC('Name',Value,...) accepts:
%     'OutDir'   destination folder                        ('figures')
%     'T'        final time                                        (6)
%     'sigma'    width of the Gaussian pulse                    (0.11)
%     'x0','y0'  centre of the pulse                          (0.3, 0)
%     'Animate'  write the animation                            (true)
%     'Npix'     animation frame size in pixels                  (460)
%     'Quick'    coarser and faster; for a first look           (false)
%     'Figures'  vector of figure numbers to (re)build          (1:7)
%
%   The solution struct is returned so the session can keep exploring:
%
%     sol = demo_wave_disc;
%     plot_solution(sol, 2.7, 'Style', 'surface');
%     sol.eval(0.1, 0.2, [0 1 2])
%
%   See also WAVE_SOLVER_DISC, PLOT_SOLUTION, ANIMATE_WAVE,
%   TEST_WAVE_SOLVER_DISC, WAVE_SOLVER_POLYGON.

p = struct('OutDir', 'figures', 'T', 6, 'sigma', 0.11, 'x0', 0.30, 'y0', 0, ...
           'Animate', true, 'Npix', 460, 'Quick', false, 'Figures', 1:7);
f = fieldnames(p);
for i = 1:2:numel(varargin)
    ix = find(strcmpi(f, varargin{i}), 1);
    if isempty(ix)
        error('demo_wave_disc:unknownOption', 'Unknown option ''%s''.', varargin{i});
    end
    p.(f{ix}) = varargin{i+1};
end

if p.Quick
    nT = 61; p.Npix = 300; p.T = min(p.T, 4); p.sigma = max(p.sigma, 0.15);
else
    nT = 241;
end

if ~exist(p.OutDir, 'dir'), mkdir(p.OutDir); end
old = set_defaults();
cleanupObj = onCleanup(@() restore_defaults(old)); %#ok<NASGU>

fprintf('\n================ demo_wave_disc ================\n');
sol = wave_solver_disc('sigma', p.sigma, 'x0', p.x0, 'y0', p.y0, ...
                       'T', p.T, 'nT', nT);

if any(p.Figures == 1), fig_evolution(sol, p);      end
if any(p.Figures == 2), fig_surfaces(sol, p);       end
if any(p.Figures == 3), fig_eigenmodes(sol, p);     end
if any(p.Figures == 4), fig_spectrum(sol, p);       end
if any(p.Figures == 5), fig_verification(sol, p);   end
if any(p.Figures == 6), fig_spacetime(sol, p);      end
if any(p.Figures == 7), fig_bc_compare(sol, p);     end

if p.Animate
    fprintf('  animation\n');
    animate_wave(sol, fullfile(p.OutDir, 'wave_disc.gif'), ...
                 'Npix', p.Npix, 'FPS', 24, 'Gamma', 0.6);
end

fprintf('================ done: see %s%s ================\n\n', p.OutDir, filesep);
end

% ======================================================================
function fig_evolution(sol, p)
%FIG_EVOLUTION  Eight snapshots of the propagating front.
fprintf('  fig1  evolution montage\n');
tv = sol.opt.T*[0 0.055 0.12 0.19 0.26 0.42 0.62 0.90];
g  = 0.6;
cl = sol.Wdisp*[-1 1];

fig = new_figure(1500, 880);
for k = 1:8
    ax = subplot_grid(fig, 2, 4, k, [0.030 0.115 0.035 0.075 0.075]);
    plot_solution(sol, tv(k), 'Axes', ax, 'CLim', cl, 'Gamma', g, ...
                  'Title', sprintf('t = %.2f', tv(k)));
end
wave_colorstrip(fig, [0.35 0.060 0.30 0.022], wave_colormap('wave'), cl, g, ...
                'displacement  w');
super_title(fig, sprintf(['Wave equation on the disc:  %s   (c = %g, R = %g, ' ...
    'Dirichlet)'], sol.icInfo.label, sol.opt.c, sol.opt.R));
save_fig(fig, fullfile(p.OutDir, 'fig1_evolution.png'));
close(fig);
end

% ======================================================================
function fig_surfaces(sol, p)
%FIG_SURFACES  The same field drawn as 3D relief.
fprintf('  fig2  3D surfaces\n');
tv = sol.opt.T*[0.12 0.30 0.62];
g  = 0.6;
cl = sol.Wdisp*[-1 1];
zl = 1.25*sol.Wdisp*[-1 1];

fig = new_figure(1500, 560);
for k = 1:3
    ax = subplot_grid(fig, 1, 3, k, [0.025 0.175 0.055 0.10 0.105]);
    plot_solution(sol, tv(k), 'Axes', ax, 'Style', 'surface', 'CLim', cl, ...
                  'Gamma', g, 'ZLim', zl, 'Lighting', true, ...
                  'Title', sprintf('t = %.2f', tv(k)));
    set(ax, 'FontSize', 9);
    zlabel(ax, 'w');
end
wave_colorstrip(fig, [0.36 0.075 0.28 0.022], wave_colormap('wave'), cl, g, ...
                'displacement  w');
super_title(fig, 'Relief view: the reflected front and its refocusing');
save_fig(fig, fullfile(p.OutDir, 'fig2_surfaces.png'));
close(fig);
end

% ======================================================================
function fig_eigenmodes(sol, p)
%FIG_EIGENMODES  Gallery of the lowest Dirichlet eigenfunctions.
%
%   These are the building blocks the solver expands in: each panel is one
%   J_m(lambda_{mn} r) cos/sin(m theta), and the label gives its angular
%   and radial indices together with the eigenvalue.
fprintf('  fig3  eigenmode gallery\n');
E = sol.E;
nsh = 20;
fig = new_figure(1400, 1180);
cmap = wave_colormap('wave');
rp = sol.rp; thc = [sol.thp, 2*pi];
for k = 1:nsh
    m = E.m(k); lam = E.lambda(k);
    if E.par(k) > 0, ang = cos(m*thc); ps = 'cos'; else, ang = sin(m*thc); ps = 'sin'; end
    Wk = besselj(m, lam*rp)*ang;
    a  = max(abs(Wk(:)));

    ax = subplot_grid(fig, 4, 5, k, [0.020 0.020 0.020 0.055 0.075]);
    hs = surf(ax, sol.X, sol.Y, zeros(size(Wk)), Wk);
    set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(ax, 2); axis(ax, 'equal'); axis(ax, 'off');
    set(ax, 'XLim', 1.04*[-1 1], 'YLim', 1.04*[-1 1], 'CLim', a*[-1 1]);
    colormap(ax, cmap);
    hold(ax, 'on');
    th = linspace(0, 2*pi, 400);
    plot3(ax, cos(th), sin(th), 1e-3*ones(size(th)), 'Color', [0.12 0.12 0.16]);
    hold(ax, 'off');
    ht = title(ax, sprintf('(%d,%d) %s   \\lambda = %.3f', m, E.n(k), ps, lam));
    set(ht, 'FontWeight', 'normal', 'FontSize', 9);
end
super_title(fig, ['Dirichlet eigenfunctions of the disc,  J_m(\lambda_{mn} r)' ...
                  ' {cos, sin}(m\theta),  ordered by \lambda']);
save_fig(fig, fullfile(p.OutDir, 'fig3_eigenmodes.png'));
close(fig);
end

% ======================================================================
function fig_spectrum(sol, p)
%FIG_SPECTRUM  Where the initial data sits in the spectrum of the disc.
fprintf('  fig4  spectral content\n');
E  = sol.E;
en = 0.5*E.nrm.*(sol.b.^2 + (sol.omega.^2).*sol.a.^2);   % modal energy
frac = en/max(sum(en), realmin);
lam = E.lambda;

fig = new_figure(1500, 520);

% (a) the modes in the (m, lambda) plane, coloured by modal energy
ax = subplot_grid(fig, 1, 3, 1, [0.055 0.135 0.085 0.10 0.135]);
cval = log10(max(frac, 1e-16));
cval = max(cval, -12);
[~, ord] = sort(cval);                      % draw the loud modes last
scatter(ax, E.m(ord), lam(ord), 16, cval(ord), 'filled');
colormap(ax, wave_colormap('ember'));
set(ax, 'CLim', [-12 max(cval)]);
cb = safe_colorbar(ax);
try, set(get(cb, 'Label'), 'String', 'log_{10} energy fraction');
catch, ylabel(cb, 'log_{10} energy fraction'); end
xlabel(ax, 'angular order  m'); ylabel(ax, 'eigenvalue  \lambda_{mn}');
title(ax, 'modal energy in the (m,\lambda) plane', 'FontWeight', 'normal');
box(ax, 'on'); grid(ax, 'on');

% (b) the radial decay, against the Gaussian envelope
ax = subplot_grid(fig, 1, 3, 2, [0.055 0.135 0.085 0.10 0.135]);
amp = abs(sol.a).*sqrt(E.nrm);
la  = log10(max(amp, 1e-18));
hsc = scatter(ax, lam, la, 12, E.m, 'filled');
set(hsc, 'HandleVisibility', 'off');      % the legend is for the envelope only
colormap(ax, wave_colormap('abyss'));
set(ax, 'YLim', [-17 ceil(max(la))]);
hold(ax, 'on');
ll = linspace(0, max(lam), 300);
env = log10(max(amp)) - sol.opt.sigma^2*ll.^2/2/log(10);
plot(ax, ll, env, '--', 'LineWidth', 1.8, 'Color', [0.75 0.13 0.13]);
hold(ax, 'off');
legend(ax, {'e^{-\sigma^2\lambda^2/2}'}, 'Location', 'northeast');
legend(ax, 'boxoff');
xlabel(ax, 'eigenvalue  \lambda');
ylabel(ax, 'log_{10} ( |a_k| \cdot ||\phi_k|| )');
title(ax, 'spectral decay of the initial data', 'FontWeight', 'normal');
box(ax, 'on'); grid(ax, 'on');
cb = safe_colorbar(ax);
try, set(get(cb, 'Label'), 'String', 'm'); catch, ylabel(cb, 'm'); end

% (c) the counting function against Weyl's two-term law
ax = subplot_grid(fig, 1, 3, 3, [0.055 0.135 0.085 0.10 0.135]);
ls = sort(lam);
stairs(ax, ls, (1:numel(ls))', 'LineWidth', 1.4, 'Color', [0.09 0.31 0.62]);
hold(ax, 'on');
R = sol.opt.R;
weyl = R^2*ll.^2/4 - R*ll/2;               % area term minus perimeter term
plot(ax, ll, max(weyl, 0), '--', 'LineWidth', 1.6, 'Color', [0.75 0.13 0.13]);
hold(ax, 'off');
xlabel(ax, '\lambda'); ylabel(ax, 'N(\lambda)');
title(ax, 'counting function vs Weyl asymptotics', 'FontWeight', 'normal');
legend(ax, {'computed basis', 'R^2\lambda^2/4 - R\lambda/2'}, 'Location', 'northwest');
legend(ax, 'boxoff'); box(ax, 'on'); grid(ax, 'on');

super_title(fig, sprintf(['Spectral content:  %d modes, \\lambda \\leq %.3g,  ' ...
    'capturing 1 - %.1e of ||w_0||'], E.K, sol.opt.LambdaMax, sol.diag.projErr));
save_fig(fig, fullfile(p.OutDir, 'fig4_spectrum.png'));
close(fig);
end

% ======================================================================
function fig_verification(sol, p)
%FIG_VERIFICATION  Convergence, conservation and accuracy of the solver.
fprintf('  fig5  verification\n');
fig = new_figure(1400, 960);
col = [0.09 0.31 0.62; 0.75 0.13 0.13; 0.10 0.53 0.35; 0.55 0.25 0.60];

% (a) spectral convergence of the projection
ax = subplot_grid(fig, 2, 2, 1, [0.060 0.075 0.085 0.115 0.085]);
lamList = [10 15 20 25 30 35 40];
pe = zeros(size(lamList)); bc = 0;
for i = 1:numel(lamList)
    si = wave_solver_disc('sigma', 0.15, 'x0', 0.15, 'LambdaMax', lamList(i), ...
                          'nT', 2, 'Verbose', false, 'PrecomputeField', false);
    pe(i) = si.diag.projErr; bc = si.diag.bcResidual;
end
semilogy(ax, lamList, pe, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6, ...
         'Color', col(1,:), 'MarkerFaceColor', col(1,:));
hold(ax, 'on');
plot(ax, lamList([1 end]), bc*[1 1], '--', 'LineWidth', 1.5, 'Color', col(2,:));
hold(ax, 'off');
xlabel(ax, 'spectral cutoff  \lambda_{max}');
ylabel(ax, 'relative L^2 projection error');
title(ax, '(a) spectral convergence of the expansion', 'FontWeight', 'normal');
legend(ax, {'projection error', 'boundary-mismatch floor'}, 'Location', 'northeast');
legend(ax, 'boxoff'); grid(ax, 'on'); box(ax, 'on');

% (b) energy conservation, closed form against ODE45
ax = subplot_grid(fig, 2, 2, 2, [0.060 0.075 0.085 0.115 0.085]);
args = {'IC', 'mode', 'modeM', 2, 'modeN', 3, 'T', 20, 'nT', 201, ...
        'Verbose', false, 'PrecomputeField', false};
sa = wave_solver_disc(args{:});
so = wave_solver_disc(args{:}, 'TimeScheme', 'ode45');
ea = modal_energy_of(sa, sa.alpha, sa.alphaDot);
eo = modal_energy_of(so, so.alpha, so.alphaDot);
semilogy(ax, sa.t, max(abs(ea/ea(1) - 1), 1e-18), 'LineWidth', 1.6, 'Color', col(1,:));
hold(ax, 'on');
semilogy(ax, so.t, max(abs(eo/eo(1) - 1), 1e-18), 'LineWidth', 1.6, 'Color', col(2,:));
hold(ax, 'off');
xlabel(ax, 't'); ylabel(ax, '|E(t)/E(0) - 1|');
title(ax, '(b) energy drift of the two time integrators', 'FontWeight', 'normal');
legend(ax, {'closed form', 'ode45 (RelTol 10^{-9})'}, 'Location', 'east');
legend(ax, 'boxoff'); grid(ax, 'on'); box(ax, 'on');

% (c) a single eigenmode against its closed form
ax = subplot_grid(fig, 2, 2, 3, [0.060 0.075 0.085 0.115 0.085]);
sm = wave_solver_disc('IC', 'mode', 'modeM', 3, 'modeN', 4, 'c', 1.3, 'T', 12, ...
                      'nT', 121, 'Verbose', false, 'PrecomputeField', false);
[xg, yg] = meshgrid(linspace(-0.92, 0.92, 21));
in = hypot(xg, yg) <= 0.92;
Wn = sm.eval(xg(in), yg(in), sm.t);
We = sm.icInfo.exact(xg(in), yg(in), sm.t);
semilogy(ax, sm.t, max(max(abs(Wn - We), [], 1), 1e-18), 'LineWidth', 1.6, ...
         'Color', col(3,:));
xlabel(ax, 't'); ylabel(ax, 'max |w_{num} - w_{exact}|');
title(ax, '(c) error against the exact single-mode solution', 'FontWeight', 'normal');
grid(ax, 'on'); box(ax, 'on');

% (d) how exactly the boundary conditions are enforced
ax = subplot_grid(fig, 2, 2, 4, [0.060 0.075 0.085 0.115 0.085]);
th = linspace(0, 2*pi, 128)';
R  = sol.opt.R;
bd = max(abs(sol.eval(R*cos(th), R*sin(th), sol.t)), [], 1);
sn = wave_solver_disc('BC', 'neumann', 'sigma', sol.opt.sigma, 'x0', sol.opt.x0, ...
                      'T', sol.opt.T, 'nT', 61, 'Verbose', false, ...
                      'PrecomputeField', false);
h  = 1e-4;
dn = (3*sn.eval(R*cos(th), R*sin(th), sn.t) ...
      - 4*sn.eval((R-h)*cos(th), (R-h)*sin(th), sn.t) ...
      + sn.eval((R-2*h)*cos(th), (R-2*h)*sin(th), sn.t))/(2*h);
semilogy(ax, sol.t, max(bd, 1e-18), 'LineWidth', 1.6, 'Color', col(1,:));
hold(ax, 'on');
semilogy(ax, sn.t, max(max(abs(dn), [], 1), 1e-18), 'LineWidth', 1.6, 'Color', col(4,:));
hold(ax, 'off');
xlabel(ax, 't'); ylabel(ax, 'residual on r = R');
title(ax, '(d) boundary conditions hold to machine precision', 'FontWeight', 'normal');
legend(ax, {'Dirichlet: max_\theta |w|', ...
            'Neumann: max_\theta |\partial_r w| (finite-difference floor)'}, ...
       'Location', 'east');
legend(ax, 'boxoff'); grid(ax, 'on'); box(ax, 'on');

super_title(fig, 'Verification of the spectral disc solver');
save_fig(fig, fullfile(p.OutDir, 'fig5_verification.png'));
close(fig);
end

% ======================================================================
function fig_spacetime(sol, p)
%FIG_SPACETIME  Space-time diagrams, and the characteristics they follow.
%
%   The left panel cuts along the diameter through the pulse; the right
%   panel uses a pulse released at the centre, whose front is a single
%   radial characteristic bouncing between r = 0 and r = R.  Overlaying
%   that characteristic shows the spectral solution reproducing finite
%   propagation speed, even though every basis function is global.
fprintf('  fig6  space-time diagrams\n');
R = sol.opt.R; c = sol.opt.c;
g = 0.55;

fig = new_figure(1400, 660);

% (a) cut along the diameter through the pulse centre
ax = subplot_grid(fig, 1, 2, 1, [0.060 0.165 0.090 0.10 0.110]);
xs  = linspace(-R, R, 421)';
tst = linspace(0, sol.opt.T, 401);
Wd  = sol.eval(xs, zeros(size(xs)), tst);
A   = sol.Wdisp;
Cd  = sign(Wd).*(min(abs(Wd)/A, 1).^g);
imagesc(ax, xs, tst, Cd');
colormap(ax, wave_colormap('wave'));
set(ax, 'YDir', 'normal', 'CLim', [-1 1], 'XLim', [-R R], 'YLim', [0 sol.opt.T]);
xlabel(ax, 'x   (along y = 0)'); ylabel(ax, 't');
title(ax, sprintf('(a) w(x,0,t),  pulse at x_0 = %.2g', sol.opt.x0), ...
      'FontWeight', 'normal');
box(ax, 'on');

% (b) a centred pulse: one radial characteristic, bouncing
ax = subplot_grid(fig, 1, 2, 2, [0.060 0.165 0.090 0.10 0.110]);
sc = wave_solver_disc('x0', 0, 'y0', 0, 'sigma', sol.opt.sigma, 'T', sol.opt.T, ...
                      'tOut', tst, 'Verbose', false, 'PrecomputeField', false);
rs = linspace(0, R, 281)';
Wr = sc.eval(rs, zeros(size(rs)), sc.t);
Ac = max(abs(Wr(:)))*0.35;
Cr = sign(Wr).*(min(abs(Wr)/Ac, 1).^g);
imagesc(ax, rs, sc.t, Cr');
colormap(ax, wave_colormap('wave'));
set(ax, 'YDir', 'normal', 'CLim', [-1 1], 'XLim', [0 R], 'YLim', [0 sol.opt.T]);
hold(ax, 'on');
% the front leaves r = 0 at speed c and folds back at r = R
tt = linspace(0, sol.opt.T, 2000);
rc = R*(1 - abs(mod(c*tt/R, 2) - 1));
plot(ax, rc, tt, '--', 'LineWidth', 1.5, 'Color', [0.10 0.10 0.12]);
hold(ax, 'off');
xlabel(ax, 'r'); ylabel(ax, 't');
title(ax, '(b) w(r,t), centred pulse, with the characteristic r = ct', ...
      'FontWeight', 'normal');
box(ax, 'on');

wave_colorstrip(fig, [0.37 0.062 0.26 0.022], wave_colormap('wave'), [-1 1], 1, ...
                'w / w_{sat},  \gamma-stretched');
super_title(fig, 'Space-time structure: reflection at r = R and refocusing');
save_fig(fig, fullfile(p.OutDir, 'fig6_spacetime.png'));
close(fig);
end

% ======================================================================
function fig_bc_compare(sol, p)
%FIG_BC_COMPARE  The boundary-condition toggle, side by side.
%
%   Same initial pulse, same times.  Under Dirichlet conditions the front
%   comes back from the rim with its sign reversed; under Neumann it comes
%   back with the sign preserved.  The contrast is clearest just after the
%   first contact with the boundary.
fprintf('  fig7  Dirichlet vs Neumann\n');
sn = wave_solver_disc('BC', 'neumann', 'sigma', sol.opt.sigma, 'x0', sol.opt.x0, ...
                      'y0', sol.opt.y0, 'T', sol.opt.T, 'nT', numel(sol.t), ...
                      'Verbose', false);
tv = sol.opt.T*[0.10 0.16 0.26];
g  = 0.6;
cl = max(sol.Wdisp, sn.Wdisp)*[-1 1];

fig = new_figure(1250, 980);
lab = {'Dirichlet   w = 0', 'Neumann   \partial w/\partial n = 0'};
for row = 1:2
    if row == 1, s = sol; else, s = sn; end
    for k = 1:3
        ax = subplot_grid(fig, 2, 3, (row-1)*3 + k, [0.055 0.115 0.030 0.075 0.075]);
        plot_solution(s, tv(k), 'Axes', ax, 'CLim', cl, 'Gamma', g, ...
                      'Title', sprintf('t = %.2f', tv(k)));
        if k == 1
            text('Parent', ax, 'Units', 'normalized', 'Position', [-0.10 0.5], ...
                 'String', lab{row}, 'Rotation', 90, 'FontSize', 12, ...
                 'HorizontalAlignment', 'center');
        end
    end
end
wave_colorstrip(fig, [0.35 0.060 0.30 0.022], wave_colormap('wave'), cl, g, ...
                'displacement  w');
super_title(fig, 'The boundary-condition toggle: reflection with and without sign reversal');
save_fig(fig, fullfile(p.OutDir, 'fig7_dirichlet_vs_neumann.png'));
close(fig);
end

% ======================================================================
%  small shared helpers
% ======================================================================
function fig = new_figure(w, h)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 w h]);
set(fig, 'PaperPositionMode', 'auto');
end

% ----------------------------------------------------------------------
function ax = subplot_grid(fig, nr, nc, k, pad)
%SUBPLOT_GRID  SUBPLOT with explicit margins, so panels can be packed
%   tightly, panel titles clear the figure title, and a shared colour key
%   can sit underneath.  PAD is [left bottom hgap vgap top] in normalised
%   figure units; the left margin is mirrored on the right.
if nargin < 5 || isempty(pad), pad = [0.05 0.08 0.05 0.10 0.08]; end
L = pad(1); B = pad(2); HG = pad(3); VG = pad(4); TP = pad(5);
Wt = (1 - 2*L - (nc-1)*HG)/nc;
Ht = (1 - B - TP - (nr-1)*VG)/nr;
r  = ceil(k/nc); cc = k - (r-1)*nc;
x  = L + (cc-1)*(Wt + HG);
y  = 1 - TP - r*Ht - (r-1)*VG;
ax = axes('Parent', fig, 'Position', [x y Wt Ht]);
end

% ----------------------------------------------------------------------
function super_title(fig, str)
%SUPER_TITLE  Title across the whole figure, without needing SGTITLE.
ax = axes('Parent', fig, 'Position', [0 0.955 1 0.04], 'Visible', 'off');
text('Parent', ax, 'Units', 'normalized', 'Position', [0.5 0.3], ...
     'String', str, 'HorizontalAlignment', 'center', 'FontSize', 13);
end

% ----------------------------------------------------------------------
function E = modal_energy_of(s, alpha, adot)
E = 0.5*sum((s.E.nrm*ones(1, size(alpha, 2))).*(adot.^2 + (s.omega.^2).*alpha.^2), 1);
end

% ----------------------------------------------------------------------
function old = set_defaults()
old = get(0, {'DefaultAxesFontSize', 'DefaultTextFontSize', ...
              'DefaultAxesLineWidth', 'DefaultAxesColor'});
set(0, 'DefaultAxesFontSize', 10, 'DefaultTextFontSize', 10, ...
       'DefaultAxesLineWidth', 0.8, 'DefaultAxesColor', 'w');
end

function restore_defaults(old)
set(0, 'DefaultAxesFontSize', old{1}, 'DefaultTextFontSize', old{2}, ...
       'DefaultAxesLineWidth', old{3}, 'DefaultAxesColor', old{4});
end
