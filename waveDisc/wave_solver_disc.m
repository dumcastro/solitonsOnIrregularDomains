function sol = wave_solver_disc(varargin)
%WAVE_SOLVER_DISC  Spectral solver for the 2D wave equation on a disc.
%
%   SOL = WAVE_SOLVER_DISC() solves
%
%       w_tt = c^2 Lap w      in  D = {(x,y): x^2 + y^2 < R^2}
%       w    = 0              on  dD          ('dirichlet', the default)
%    or dw/dn = 0             on  dD          ('neumann')
%       w(.,0) = w0,   w_t(.,0) = v0
%
%   by separation of variables in polar coordinates.  The solution is
%   expanded in the exact Bessel eigenbasis of the disc,
%
%       w(r,theta,t) = SUM_k alpha_k(t) phi_k(r,theta),
%       phi_k = J_m(lambda_k r) * {cos(m theta), sin(m theta)},
%
%   which turns the PDE into K decoupled harmonic oscillators
%
%       alpha_k'' + omega_k^2 alpha_k = 0,   omega_k = c*lambda_k,
%
%   integrated either in closed form (default) or with ODE45.
%
%   SOL = WAVE_SOLVER_DISC('Name',Value,...) accepts:
%
%   Geometry and physics
%     'R'          disc radius                                      (1)
%     'c'          wave speed                                       (1)
%     'BC'         'dirichlet' | 'neumann'                 ('dirichlet')
%
%   Initial data (see INITIAL_CONDITION for the full list)
%     'IC'         'gaussian' | 'ring' | 'twobumps' | 'mode' | 'custom'
%                                                          ('gaussian')
%     'amplitude'  peak amplitude of the profile                     (1)
%     'x0','y0'    Gaussian centre                              (0.3, 0)
%     'sigma'      Gaussian width                                (0.11)
%     'r0'         radius of the 'ring' profile                   (0.5)
%     'modeM','modeN','modePar'   eigenmode indices for IC 'mode'
%                                                            (2, 3, +1)
%     'w0fun'      handle @(x,y) for IC 'custom'                    ([])
%     'v0fun'      handle @(x,y) for the initial velocity           ([])
%
%   Time
%     'T'          final time                                       (4)
%     'nT'         number of output times in [0,T]                (161)
%     'tOut'       explicit vector of output times, overrides T/nT ([])
%     'TimeScheme' 'analytic' | 'ode45'                    ('analytic')
%     'RelTol','AbsTol'  ODE45 tolerances                (1e-9, 1e-11)
%
%   Discretisation
%     'LambdaMax'  spectral cutoff: keep all lambda_k <= LambdaMax
%                  (default: chosen from sigma so the data is resolved
%                   to ~1e-14, capped by MaxModes)
%     'MaxModes'   cap on the number of basis functions           (4000)
%     'Nquad'      Gauss-Legendre nodes in r for the projection   (auto)
%     'Ntheta'     angular quadrature points (FFT)                (auto)
%     'PlotNr','PlotNtheta'  polar grid used for the stored field
%                                                            (160, 384)
%     'PrecomputeField'  store w on the plot grid at all output
%                        times (needed by PLOT_SOLUTION)          (true)
%     'Verbose'    print a short report                           (true)
%
%   The returned struct SOL contains
%     .t            1-by-nt output times
%     .E            the eigenbasis (see DISC_EIGENMODES)
%     .a, .b        modal coefficients of w0 and v0
%     .omega        modal frequencies c*lambda
%     .alpha        K-by-nt modal amplitudes at SOL.t
%     .X, .Y, .W    polar plot grid (Nr-by-Nth) and field (Nr-by-Nth-by-nt)
%     .Wmax         peak |w| over the whole run
%     .Wdisp        robust display amplitude: the 99.5th percentile of |w|
%                   over all output times after t = 0.  The initial pulse is
%                   far taller than anything the spread-out wave reaches, so
%                   scaling colours to .Wmax washes the animation out;
%                   PLOT_SOLUTION therefore saturates at .Wdisp by default
%     .rp, .thp     the underlying polar grid vectors
%     .eval(x,y,t)  solution at arbitrary points and times
%     .energy(t)    total energy 0.5*INT (w_t^2 + c^2|grad w|^2)
%     .diag         resolution diagnostics (see below)
%     .opt          the resolved option struct
%
%   Note SOL.alpha and SOL.W follow the chosen TimeScheme, whereas
%   SOL.eval, SOL.field and SOL.energy always use the closed-form modal
%   propagator so that they can be queried at arbitrary times.
%
%   SOL.diag reports .captured (fraction of the L2 energy of w0 represented
%   in the truncated basis), .bcResidual (how badly w0 violates the boundary
%   condition, relative to its peak) and .projErr (L2 error of the projected
%   initial data on the quadrature grid).  All three should be tiny; the
%   solver warns when they are not.
%
%   Example
%     sol = wave_solver_disc('sigma',0.1,'x0',0.4,'T',6);
%     plot_solution(sol, 1.5, 'Style','surface');
%
%   See also INITIAL_CONDITION, PLOT_SOLUTION, DISC_EIGENMODES,
%   DEMO_WAVE_DISC, WAVE_SOLVER_POLYGON.

% ---------------------------------------------------------------- options
d = struct( ...
    'R', 1, 'c', 1, 'BC', 'dirichlet', ...
    'IC', 'gaussian', 'amplitude', 1, ...
    'x0', 0.30, 'y0', 0, 'sigma', 0.11, 'r0', 0.5, ...
    'modeM', 2, 'modeN', 3, 'modePar', 1, ...
    'w0fun', [], 'v0fun', [], ...
    'T', 4, 'nT', 161, 'tOut', [], ...
    'TimeScheme', 'analytic', 'RelTol', 1e-9, 'AbsTol', 1e-11, ...
    'LambdaMax', [], 'MaxModes', 4000, 'Nquad', [], 'Ntheta', [], ...
    'PlotNr', 160, 'PlotNtheta', 384, 'PrecomputeField', true, ...
    'Verbose', true);

opt = parse_opts(d, varargin);

if opt.R <= 0, error('wave_solver_disc:badR', 'R must be positive.'); end
if opt.c <= 0, error('wave_solver_disc:badC', 'c must be positive.'); end
if ~strncmpi(opt.BC, 'd', 1) && ~strncmpi(opt.BC, 'n', 1)
    error('wave_solver_disc:badBC', 'BC must be ''dirichlet'' or ''neumann''.');
end
if strcmpi(opt.IC, 'gaussian') && hypot(opt.x0, opt.y0) >= opt.R
    error('wave_solver_disc:centreOutside', ...
        'The Gaussian centre (%.3g, %.3g) lies outside the disc of radius %.3g.', ...
        opt.x0, opt.y0, opt.R);
end

if isempty(opt.tOut)
    t = linspace(0, opt.T, opt.nT);
else
    t = opt.tOut(:)';
    opt.T = max(t);
end

% ------------------------------------------------------- initial data
[w0fun, v0fun, icInfo] = initial_condition(opt.IC, opt);

% ------------------------------------------------------- spectral cutoff
if isempty(opt.LambdaMax)
    opt.LambdaMax = auto_lambda_max(opt, icInfo);
end

E = disc_eigenmodes(opt.R, opt.BC, opt.LambdaMax, opt.MaxModes);
K = E.K;
omega = opt.c*E.lambda;

% ------------------------------------------------- quadrature for the
%                                                    L2 projection
if isempty(opt.Nquad)
    Nr = ceil(1.6*opt.LambdaMax*opt.R) + 80;
else
    Nr = opt.Nquad;
end
if isempty(opt.Ntheta)
    Nth = 2^nextpow2(4*E.Mmax + 32);
else
    Nth = opt.Ntheta;
end
if Nth <= 2*E.Mmax
    error('wave_solver_disc:thetaAliasing', ...
        'Ntheta = %d aliases angular orders up to m = %d; use Ntheta > %d.', ...
        Nth, E.Mmax, 2*E.Mmax);
end

[rq, wq] = gauss_legendre(Nr, 0, opt.R);
thq = 2*pi*(0:Nth-1)/Nth;
Xq = rq*cos(thq);            % Nr-by-Nth, outer products
Yq = rq*sin(thq);

F0 = w0fun(Xq, Yq);
V0 = v0fun(Xq, Yq);
if isscalar(V0), V0 = V0*ones(size(Xq)); end
if ~isequal(size(F0), size(Xq))
    error('wave_solver_disc:icShape', ...
        'The initial displacement handle must return an array the size of its input.');
end

% Angular transform: Fh(:,m+1) = INT_0^{2pi} f(r,th) exp(-i m th) dth.
Fh = fft(F0, [], 2)*(2*pi/Nth);
Vh = fft(V0, [], 2)*(2*pi/Nth);

% ------------------------------------------------------- mode bookkeeping
mList = unique(E.m(:))';
blk = struct('m', {}, 'idxC', {}, 'idxS', {}, 'lam', {});
for m = mList
    b.m    = m;
    b.idxC = find(E.m == m & E.par > 0);
    b.idxS = find(E.m == m & E.par < 0);
    b.lam  = E.lambda(b.idxC);
    blk(end+1) = b; %#ok<AGROW>
end

% ------------------------------------------------------- L2 projection
a = zeros(K, 1);   % coefficients of w0
b_ = zeros(K, 1);  % coefficients of v0
rw = wq.*rq;       % radial quadrature weight including the Jacobian r
for q = 1:numel(blk)
    m  = blk(q).m;
    B  = besselj(m, rq*blk(q).lam');            % Nr-by-nm
    gC = real(Fh(:, m+1));   hC = real(Vh(:, m+1));
    a (blk(q).idxC) = (B'*(rw.*gC))./E.nrm(blk(q).idxC);
    b_(blk(q).idxC) = (B'*(rw.*hC))./E.nrm(blk(q).idxC);
    if ~isempty(blk(q).idxS)
        gS = -imag(Fh(:, m+1));  hS = -imag(Vh(:, m+1));
        a (blk(q).idxS) = (B'*(rw.*gS))./E.nrm(blk(q).idxS);
        b_(blk(q).idxS) = (B'*(rw.*hS))./E.nrm(blk(q).idxS);
    end
end

% ------------------------------------------------------- diagnostics
dth      = 2*pi/Nth;
normW0sq = sum(rw'*(F0.^2))*dth;
captured = sum(E.nrm.*a.^2);
dg = struct();
dg.normW0   = sqrt(normW0sq);
dg.captured = captured/max(normW0sq, realmin);
dg.projErr  = sqrt(max(normW0sq - captured, 0)/max(normW0sq, realmin));
dg.nModes   = K;
dg.Mmax     = E.Mmax;
dg.Nquad    = Nr;
dg.Ntheta   = Nth;

thb = linspace(0, 2*pi, 721);
if strncmpi(opt.BC, 'd', 1)
    bres = max(abs(w0fun(opt.R*cos(thb), opt.R*sin(thb))));
    dg.bcType = 'w on dD';
else
    h    = 1e-5*opt.R;
    wIn  = w0fun((opt.R - h)*cos(thb), (opt.R - h)*sin(thb));
    wIn2 = w0fun((opt.R - 2*h)*cos(thb), (opt.R - 2*h)*sin(thb));
    wOn  = w0fun(opt.R*cos(thb), opt.R*sin(thb));
    bres = max(abs((3*wOn - 4*wIn + wIn2)/(2*h)))*opt.R;   % R*dw/dn
    dg.bcType = 'R*dw/dn on dD';
end
peak = max(abs(F0(:)));
dg.bcResidual = bres/max(peak, realmin);

if dg.projErr > 1e-6
    warning('wave_solver_disc:underResolved', ...
        ['The truncated basis captures only 1 - %.2e of the initial data.\n' ...
         'Raise LambdaMax (currently %.4g) and MaxModes, or widen sigma.'], ...
        dg.projErr, opt.LambdaMax);
end
if dg.bcResidual > 1e-6
    warning('wave_solver_disc:boundaryMismatch', ...
        ['The initial data violates the %s boundary condition by %.2e of its peak\n' ...
         '(%s). The eigenfunction expansion then converges slowly (Gibbs); move\n' ...
         'the pulse away from r = R or narrow sigma.'], ...
        E.bc, dg.bcResidual, dg.bcType);
end

% ------------------------------------------------------- time evolution
switch lower(opt.TimeScheme)
    case 'analytic'
        alpha = modal_analytic(a, b_, omega, t);
        adot  = modal_analytic_dot(a, b_, omega, t);
    case 'ode45'
        [alpha, adot] = modal_ode45(a, b_, omega, t, opt);
    otherwise
        error('wave_solver_disc:badScheme', ...
            'TimeScheme must be ''analytic'' or ''ode45''.');
end

% ------------------------------------------------------- output struct
sol = struct();
sol.opt    = opt;
sol.E      = E;
sol.t      = t;
sol.a      = a;
sol.b      = b_;
sol.omega  = omega;
sol.alpha  = alpha;
sol.alphaDot = adot;
sol.diag   = dg;
sol.icInfo = icInfo;
sol.w0fun  = w0fun;
sol.v0fun  = v0fun;
sol.blocks = blk;

sol.eval   = @(x, y, tt) eval_points(blk, E, a, b_, omega, x, y, tt);
sol.energy = @(tt) modal_energy(E.nrm, a, b_, omega, tt);
sol.E0     = sol.energy(0);

% Polar plotting grid: r includes 0 and R, theta is uniform and closed.
Nrp  = opt.PlotNr;
Nthp = max(opt.PlotNtheta, 2*E.Mmax + 2);
rp   = linspace(0, opt.R, Nrp)';
thp  = 2*pi*(0:Nthp-1)/Nthp;
sol.rp  = rp;
sol.thp = thp;
sol.X   = rp*cos([thp, 2*pi]);          % closed in theta for seamless plots
sol.Y   = rp*sin([thp, 2*pi]);

if opt.PrecomputeField
    sol.W = synth_polar(blk, rp, Nthp, alpha);
    sol.W = cat(2, sol.W, sol.W(:, 1, :));    % close the seam
    sol.Wmax = max(abs(sol.W(:)));
    sol.Wdisp = display_scale(sol.W);
else
    sol.W = [];
    sol.Wmax = max(abs(F0(:)));
    sol.Wdisp = sol.Wmax;
end
sol.field = @(tt) close_seam(synth_polar(blk, rp, Nthp, ...
                    modal_analytic(a, b_, omega, tt)));

% ------------------------------------------------------- report
if opt.Verbose
    fprintf('\n  wave_solver_disc  --  %s BC, R = %g, c = %g\n', E.bc, opt.R, opt.c);
    fprintf('  initial data : %s\n', icInfo.label);
    fprintf('  basis        : %d modes, lambda <= %.4g (m <= %d), omega_max = %.4g\n', ...
        K, opt.LambdaMax, E.Mmax, max(omega));
    fprintf('  quadrature   : %d Gauss-Legendre nodes in r, %d in theta\n', Nr, Nth);
    fprintf('  time         : %s on [0, %g], %d output times\n', ...
        lower(opt.TimeScheme), opt.T, numel(t));
    fprintf('  projection   : captures 1 - %.2e of ||w0||_L2\n', dg.projErr);
    fprintf('  BC residual  : %.2e of peak (%s)\n', dg.bcResidual, dg.bcType);
    fprintf('  energy       : E(0) = %.10g, drift over [0,T] = %.2e\n\n', ...
        sol.E0, max(abs(sol.energy(t) - sol.E0))/max(sol.E0, realmin));
end

end % wave_solver_disc

% ======================================================================
function alpha = modal_analytic(a, b, omega, t)
%MODAL_ANALYTIC  Exact solution of alpha'' + omega^2 alpha = 0.
t  = t(:)';
ct = cos(omega*t);
alpha = (a*ones(1, numel(t))).*ct;
nz = omega > 0;
if any(nz)
    alpha(nz, :) = alpha(nz, :) + (b(nz)./omega(nz)).*sin(omega(nz)*t);
end
if any(~nz)                              % Neumann constant mode: a + b*t
    alpha(~nz, :) = a(~nz)*ones(1, numel(t)) + b(~nz)*t;
end
end

% ----------------------------------------------------------------------
function adot = modal_analytic_dot(a, b, omega, t)
t = t(:)';
adot = zeros(numel(a), numel(t));
nz = omega > 0;
adot(nz, :) = -(a(nz).*omega(nz)).*sin(omega(nz)*t) + b(nz).*cos(omega(nz)*t);
adot(~nz, :) = b(~nz)*ones(1, numel(t));
end

% ----------------------------------------------------------------------
function [alpha, adot] = modal_ode45(a, b, omega, t, opt)
%MODAL_ODE45  Same modal system, integrated with MATLAB's ODE45.
%
%   Kept because the modal system y' = [ydot; -omega.^2 .* y] is the
%   natural place to swap in a different time integrator (a forcing term,
%   damping, or the variable-coefficient polygon problem of Phase 2), where
%   the closed form is no longer available.
w2 = omega.^2;
rhs = @(tt, y) [y(numel(a)+1:end); -w2.*y(1:numel(a))];
odeopt = odeset('RelTol', opt.RelTol, 'AbsTol', opt.AbsTol);
tt = t(:);
flip = false;
if numel(tt) < 2                          % ODE45 needs a real interval
    tt = [tt; tt + 1]; flip = true;
end
[~, Y] = ode45(rhs, tt, [a(:); b(:)], odeopt);
if flip, Y = Y(1, :); end
alpha = Y(:, 1:numel(a))';
adot  = Y(:, numel(a)+1:end)';
end

% ----------------------------------------------------------------------
function W = synth_polar(blk, rp, Nth, alpha)
%SYNTH_POLAR  Field on a uniform polar grid by inverse FFT in theta.
%
%   w(r,th) = Re SUM_m ( C_m(r) - i S_m(r) ) exp(i m th), and with
%   th_j = 2 pi j / Nth that sum is exactly Nth*IFFT of the coefficient
%   array along the angular dimension.
nt  = size(alpha, 2);
Nrp = numel(rp);
Bc  = cell(1, numel(blk));
for q = 1:numel(blk)
    Bc{q} = besselj(blk(q).m, rp(:)*blk(q).lam');
end
W = zeros(Nrp, Nth, nt);
for it = 1:nt
    V = zeros(Nrp, Nth);
    for q = 1:numel(blk)
        m = blk(q).m;
        z = alpha(blk(q).idxC, it);
        if ~isempty(blk(q).idxS)
            z = z - 1i*alpha(blk(q).idxS, it);
        end
        V(:, m+1) = Bc{q}*z;
    end
    W(:, :, it) = real(Nth*ifft(V, [], 2));
end
end

% ----------------------------------------------------------------------
function W = close_seam(W)
W = cat(2, W, W(:, 1, :));
end

% ----------------------------------------------------------------------
function W = eval_points(blk, E, a, b, omega, x, y, t) %#ok<INUSL>
%EVAL_POINTS  Solution at arbitrary Cartesian points and times.
sz = size(x);
x = x(:); y = y(:);
r  = hypot(x, y);
th = atan2(y, x);
alpha = modal_analytic(a, b, omega, t);
nt = size(alpha, 2);
W  = zeros(numel(x), nt);
for q = 1:numel(blk)
    m = blk(q).m;
    B = besselj(m, r*blk(q).lam');
    W = W + (B*alpha(blk(q).idxC, :)).*cos(m*th);
    if ~isempty(blk(q).idxS)
        W = W + (B*alpha(blk(q).idxS, :)).*sin(m*th);
    end
end
if nt == 1
    W = reshape(W, sz);
end
end

% ----------------------------------------------------------------------
function Ev = modal_energy(nrm, a, b, omega, t)
%MODAL_ENERGY  0.5*INT (w_t^2 + c^2 |grad w|^2) = 0.5 SUM nrm (ad^2 + w^2 al^2).
alpha = modal_analytic(a, b, omega, t);
adot  = modal_analytic_dot(a, b, omega, t);
Ev = 0.5*sum((nrm*ones(1, size(alpha, 2))).*(adot.^2 + (omega.^2).*alpha.^2), 1);
end

% ----------------------------------------------------------------------
function A = display_scale(W)
%DISPLAY_SCALE  Robust symmetric colour amplitude for the whole run.
%
%   The initial pulse is a tall narrow spike; the wave that spreads from it
%   is several times smaller, so the peak is a poor colour scale.  Take the
%   99.5th percentile of |w| over every output time except t = 0, which
%   saturates only the brief refocusing events.  Uses SORT rather than
%   QUANTILE so that no toolbox is needed.
if size(W, 3) > 1
    V = abs(W(:, :, 2:end));
else
    V = abs(W);
end
V = sort(V(:));
if isempty(V) || V(end) == 0
    A = 1; return
end
A = V(max(1, ceil(0.995*numel(V))));
A = min(max(A, 0.02*V(end)), V(end));
end

% ----------------------------------------------------------------------
function lam = auto_lambda_max(opt, icInfo)
%AUTO_LAMBDA_MAX  Cutoff that resolves the initial data to ~1e-14.
%
%   The radial transform of a Gaussian of width sigma decays like
%   exp(-sigma^2 lambda^2/2), so lambda = 8/sigma leaves a relative tail of
%   about exp(-32) ~ 1e-14.  Weyl's law N(lambda) ~ R^2 lambda^2/4 then
%   converts the MaxModes budget into a hard ceiling on lambda.
if strcmpi(opt.IC, 'mode') && isfield(icInfo, 'lambda')
    lam = max(1.5*icInfo.lambda, 20/opt.R);
elseif strcmpi(opt.IC, 'custom')
    lam = 2*sqrt(opt.MaxModes)/opt.R;
else
    lam = 8/max(opt.sigma, eps);
end
cap = 2*sqrt(opt.MaxModes)/opt.R;
lam = min(max(lam, 20/opt.R), cap);
end

% ----------------------------------------------------------------------
function opt = parse_opts(opt, args)
%PARSE_OPTS  Name/value parsing, also accepting a single option struct.
if numel(args) == 1 && isstruct(args{1})
    f = fieldnames(args{1});
    for i = 1:numel(f)
        opt = assign(opt, f{i}, args{1}.(f{i}));
    end
    return
end
if mod(numel(args), 2) ~= 0
    error('wave_solver_disc:badArgs', ...
        'Options must be given as Name/Value pairs (or one option struct).');
end
for i = 1:2:numel(args)
    if ~ischar(args{i})
        error('wave_solver_disc:badArgs', 'Option %d is not a name.', (i+1)/2);
    end
    opt = assign(opt, args{i}, args{i+1});
end
end

function opt = assign(opt, name, value)
f  = fieldnames(opt);
ix = find(strcmpi(f, name), 1);
if isempty(ix)
    error('wave_solver_disc:unknownOption', ...
        'Unknown option ''%s''. Type ''help wave_solver_disc'' for the list.', name);
end
opt.(f{ix}) = value;
end
