function sol = wave_solver_polygon(varargin)
%WAVE_SOLVER_POLYGON  Phase 2 skeleton: the wave equation on a polygon,
%   solved in the disc basis through a conformal map.
%
%   SOL = WAVE_SOLVER_POLYGON('Vertices',V,...) solves
%
%       w_tt = c^2 Lap_z w   in a polygon Omega,   w = 0 (or dw/dn = 0)
%
%   by pulling Omega back to the unit disc with a Schwarz-Christoffel map
%   and reusing the Bessel machinery of WAVE_SOLVER_DISC.
%
%   THE POINT TO GET RIGHT.  A conformal map z = f(zeta) does not turn the
%   wave equation on Omega into the wave equation on the disc.  Under
%   z = f(zeta) the Laplacian picks up the conformal factor,
%
%       Lap_z = |f'(zeta)|^{-2} Lap_zeta,
%
%   so the pulled-back problem is the variable-coefficient equation
%
%       rho(zeta) w_tt = c^2 Lap_zeta w,      rho = |f'(zeta)|^2,
%
%   on the unit disc.  Only the Laplace and Helmholtz problems are
%   conformally clean; the wave equation is not.  What survives is that
%   separation of variables still works, with a WEIGHTED eigenproblem
%
%       -Lap_zeta phi = mu^2 rho(zeta) phi   in D,   phi = 0 on dD,
%
%   after which w = SUM_k c_k(t) phi_k with c_k'' + (c mu_k)^2 c_k = 0,
%   exactly as on the disc.  So the disc solver is reused not as a
%   black box but as a BASIS: this routine assembles the weighted problem
%   in the Bessel basis and diagonalises it.
%
%   Note also that rho is the Jacobian of the map, so INT_D rho u v is the
%   L2(Omega) inner product written in disc coordinates: the weighted mass
%   matrix below is exactly the Gram matrix of the pulled-back Bessel
%   functions in the physical L2, and no extra change of measure is needed
%   anywhere else.
%
%   METHOD.  With {phi_j} the Dirichlet (or Neumann) Bessel basis of the
%   disc supplied by DISC_EIGENMODES,
%
%       S_ij = INT grad phi_i . grad phi_j = delta_ij lambda_i^2 ||phi_i||^2
%       M_ij = INT rho phi_i phi_j                     (dense, quadrature)
%
%   and the generalised symmetric problem S v = mu^2 M v gives the polygon
%   eigenpairs.  S is diagonal because the Bessel functions already
%   diagonalise the Dirichlet Laplacian on the disc.
%
%   OPTIONS
%     'Vertices'   polygon vertices, complex or an N-by-2 array.  Requires
%                  Driscoll's Schwarz-Christoffel toolbox on the path.
%     'MapFun'     handle @(zeta) z, a conformal map of the unit disc, used
%                  instead of 'Vertices'.  Together with 'MapDiff' this
%                  lets the framework be exercised, and verified, with no
%                  toolbox present.
%     'MapDiff'    handle @(zeta) f'(zeta).  If omitted it is obtained by
%                  a centred difference on MapFun.
%     'c'          wave speed                                        (1)
%     'BC'         'dirichlet' | 'neumann'                 ('dirichlet')
%     'IC'         'gaussian' | 'custom' - the profile is given in the
%                  PHYSICAL plane and pulled back              ('gaussian')
%     'x0','y0','sigma','amplitude'   Gaussian parameters, in z
%     'T','nT'     time span and number of output times          (4, 121)
%     'LambdaMax'  cutoff of the disc basis.  Cost is O(K^2 Nq) to
%                  assemble M and O(K^3) to diagonalise, so this is the
%                  knob that matters.                               (14)
%     'MaxModes'   cap on the basis size                            (600)
%     'PlotNr','PlotNtheta'   grid pushed forward for plotting  (120, 256)
%     'Verbose'    print a short report                           (true)
%
%   SOL carries .mu (polygon eigenvalues), .V (disc-basis coefficients of
%   each polygon eigenfunction), .X, .Y (the plot grid pushed forward to
%   Omega), .W (the field there at SOL.t), .boundary (the image of the unit
%   circle, so PLOT_SOLUTION can outline the domain) and .t, .a, .b,
%   .alpha, .energy as in WAVE_SOLVER_DISC.
%
%   STATUS.  Working but experimental, and deliberately the simplest thing
%   that is correct.  Known limits and the natural next steps:
%     * The dense M costs O(K^2 Nq); a fast transform in theta would cut
%       that substantially, since rho is smooth.
%     * Schwarz-Christoffel maps crowd badly for elongated polygons, and
%       rho then varies over orders of magnitude, which both stiffens the
%       eigenproblem and starves the quadrature.  Reasonable aspect ratios
%       only, as the brief says.
%     * Reentrant corners make grad w singular there; the Bessel basis
%       converges algebraically rather than spectrally in that case.
%     * The Neumann path goes through the same generalised solve and keeps
%       the zero mode, but has had far less exercise than Dirichlet.
%
%   Example (no toolbox needed: an analytic map onto a rounded triangle)
%     f  = @(z) z + 0.22*z.^4;
%     df = @(z) 1 + 0.88*z.^3;
%     sol = wave_solver_polygon('MapFun', f, 'MapDiff', df, 'T', 3, ...
%                               'LambdaMax', 26, 'sigma', 0.3);
%     plot_solution(sol, 1.2);
%
%   See also WAVE_SOLVER_DISC, DISC_EIGENMODES, TEST_WAVE_SOLVER_POLYGON,
%   DEMO_WAVE_POLYGON.

d = struct('Vertices', [], 'MapFun', [], 'MapDiff', [], ...
           'c', 1, 'BC', 'dirichlet', 'IC', 'gaussian', ...
           'amplitude', 1, 'x0', 0.25, 'y0', 0, 'sigma', 0.18, ...
           'w0fun', [], 'v0fun', [], ...
           'T', 4, 'nT', 121, 'LambdaMax', 14, 'MaxModes', 600, ...
           'Nquad', [], 'Ntheta', [], 'PlotNr', 120, 'PlotNtheta', 256, ...
           'Verbose', true);
opt = parse_opts(d, varargin);

% ------------------------------------------------------- the map
[f, df, mapName] = resolve_map(opt);

% ------------------------------------------------------- disc basis
E = disc_eigenmodes(1, opt.BC, opt.LambdaMax, opt.MaxModes);
K = E.K;

% ------------------------------------------------------- quadrature
if isempty(opt.Nquad),  Nr  = ceil(3*opt.LambdaMax) + 60;        else, Nr  = opt.Nquad;  end
if isempty(opt.Ntheta), Nth = 2^nextpow2(4*E.Mmax + 32);         else, Nth = opt.Ntheta; end
[rq, wq] = gauss_legendre(Nr, 0, 1);
thq = 2*pi*(0:Nth-1)/Nth;
ZT  = rq*exp(1i*thq);                       % Nr-by-Nth points of the disc
qw  = (wq.*rq)*ones(1, Nth)*(2*pi/Nth);     % quadrature weights, with r dA

rho = abs(df(ZT)).^2;
if any(~isfinite(rho(:)))
    error('wave_solver_polygon:badMap', ...
        'The map derivative is not finite on the quadrature grid.');
end
if min(rho(:)) <= 0
    error('wave_solver_polygon:notConformal', ...
        'f'' vanishes inside the disc, so the map is not conformal there.');
end

% ------------------------------------------------------- basis on the grid
Phi = basis_matrix(E, rq, thq);              % (Nr*Nth)-by-K
wcol = qw(:);
rcol = rho(:);

% ------------------------------------------------------- S and M
S = diag(E.lambda.^2 .* E.nrm);              % exactly diagonal
M = Phi'*((wcol.*rcol)*ones(1, K).*Phi);
M = 0.5*(M + M');                            % symmetrise against roundoff

% ------------------------------------------------------- eigenproblem
[V, D] = eig(S, M);
mu2 = real(diag(D));
[mu2, ix] = sort(mu2);
V   = real(V(:, ix));
mu2 = max(mu2, 0);
mu  = sqrt(mu2);

% M-orthonormalise, so that the physical L2 norms are 1
for k = 1:K
    nk = sqrt(max(V(:, k)'*M*V(:, k), realmin));
    V(:, k) = V(:, k)/nk;
end
omega = opt.c*mu;

% ------------------------------------------------------- initial data
[w0fun, v0fun] = initial_condition(opt.IC, opt);
Z0 = f(ZT);                                  % physical positions of the nodes
F0 = w0fun(real(Z0), imag(Z0));
V0 = v0fun(real(Z0), imag(Z0));
if isscalar(V0), V0 = V0*ones(size(F0)); end

% Coefficients in the M-orthonormal polygon eigenbasis are plain physical
% L2 inner products: a_k = INT_Omega w0 psi_k = INT_D rho w0 phi_k.
rhs0 = Phi'*(wcol.*rcol.*F0(:));
rhsv = Phi'*(wcol.*rcol.*V0(:));
a = V'*rhs0;
b = V'*rhsv;

t = linspace(0, opt.T, opt.nT);
alpha = modal_analytic(a, b, omega, t);

% ------------------------------------------------------- plot grid
rp  = linspace(0, 1, opt.PlotNr)';
thp = 2*pi*(0:opt.PlotNtheta-1)/opt.PlotNtheta;
thc = [thp, 2*pi];
Zp  = rp*exp(1i*thc);
Wz  = f(Zp);
Pp  = basis_matrix(E, rp, thc);              % (Nrp*Nthc)-by-K
Cp  = Pp*(V*alpha);                          % values at the plot nodes
W   = reshape(Cp, [numel(rp), numel(thc), numel(t)]);

thb = linspace(0, 2*pi, 800);
bnd = f(exp(1i*thb));

sol = struct();
sol.opt   = opt;
sol.E     = E;
sol.map   = f;
sol.mapd  = df;
sol.mapName = mapName;
sol.mu    = mu;
sol.omega = omega;
sol.V     = V;
sol.M     = M;
sol.t     = t;
sol.a     = a;
sol.b     = b;
sol.alpha = alpha;
sol.X     = real(Wz);
sol.Y     = imag(Wz);
sol.rp    = rp;
sol.thp   = thp;
sol.W     = W;
sol.Wmax  = max(abs(W(:)));
sol.Wdisp = display_scale(W);
sol.boundary = [real(bnd(:)), imag(bnd(:))];
sol.energy   = @(tt) 0.5*sum((alpha_of(a, b, omega, tt, 'dot').^2 ...
                             + (omega.^2).*alpha_of(a, b, omega, tt, 'val').^2), 1);
sol.E0    = sol.energy(0);
sol.field = @(tt) reshape(Pp*(V*modal_analytic(a, b, omega, tt)), ...
                          [numel(rp), numel(thc), numel(tt)]);
sol.evalDisc = @(zeta, tt) basis_matrix(E, abs(zeta(:)), 0, angle(zeta(:))) ...
                           *(V*modal_analytic(a, b, omega, tt));

if opt.Verbose
    fprintf('\n  wave_solver_polygon  --  %s BC, c = %g\n', E.bc, opt.c);
    fprintf('  map          : %s\n', mapName);
    fprintf('  disc basis   : %d modes, lambda <= %.3g\n', K, opt.LambdaMax);
    fprintf('  quadrature   : %d x %d, conformal factor rho in [%.3g, %.3g]\n', ...
        Nr, Nth, min(rho(:)), max(rho(:)));
    fprintf('  spectrum     : mu_1..mu_4 = %.5f %.5f %.5f %.5f\n', mu(1:min(4,K)));
    fprintf('  energy       : E(0) = %.8g, drift over [0,T] = %.2e\n\n', ...
        sol.E0, max(abs(sol.energy(t) - sol.E0))/max(sol.E0, realmin));
end
end

% ======================================================================
function [f, df, name] = resolve_map(opt)
%RESOLVE_MAP  Either an explicit analytic map, or a Schwarz-Christoffel one.
if ~isempty(opt.MapFun)
    f = opt.MapFun;
    if ~isempty(opt.MapDiff)
        df = opt.MapDiff;
    else
        h  = 1e-6;
        df = @(z) (f(z + h) - f(z - h))/(2*h);
    end
    name = 'analytic map supplied by the caller';
    return
end

if isempty(opt.Vertices)
    error('wave_solver_polygon:noDomain', ...
        'Supply either ''Vertices'' (needs the SC toolbox) or ''MapFun''.');
end

if exist('diskmap', 'file') ~= 2 || exist('polygon', 'file') ~= 2
    error('wave_solver_polygon:noSCToolbox', ...
        ['''Vertices'' needs Driscoll''s Schwarz-Christoffel toolbox on the\n' ...
         'MATLAB path (it supplies POLYGON and DISKMAP).  Get it from\n' ...
         'https://github.com/tobydriscoll/sc-toolbox and ADDPATH it, or pass\n' ...
         'an analytic map through ''MapFun'' instead.']);
end

v = opt.Vertices;
if ~isvector(v) || isreal(v)
    v = v(:, 1) + 1i*v(:, 2);
end
P  = polygon(v(:));
mp = diskmap(P);
mp = center(mp, mean(v(:)));         % put the disc centre near the middle
f  = @(z) eval(mp, z);
df = @(z) evaldiff(mp, z);
name = sprintf('Schwarz-Christoffel disc map onto a %d-gon', numel(v));
end

% ----------------------------------------------------------------------
function Phi = basis_matrix(E, r, th, thDirect)
%BASIS_MATRIX  Values of every basis function at a set of points.
%
%   With TH a row vector the points are the tensor grid r x th, flattened
%   columnwise; with THDIRECT given, R and THDIRECT are paired pointwise.
r = r(:);
if nargin >= 4 && ~isempty(thDirect)
    th = thDirect(:);
    np = numel(r);
    Phi = zeros(np, E.K);
    for k = 1:E.K
        B = besselj(E.m(k), E.lambda(k)*r);
        if E.par(k) > 0, Phi(:, k) = B.*cos(E.m(k)*th);
        else,            Phi(:, k) = B.*sin(E.m(k)*th);
        end
    end
    return
end
th = th(:)';
np = numel(r)*numel(th);
Phi = zeros(np, E.K);
for k = 1:E.K
    B = besselj(E.m(k), E.lambda(k)*r);
    if E.par(k) > 0, A = B*cos(E.m(k)*th); else, A = B*sin(E.m(k)*th); end
    Phi(:, k) = A(:);
end
end

% ----------------------------------------------------------------------
function alpha = modal_analytic(a, b, omega, t)
t  = t(:)';
alpha = (a*ones(1, numel(t))).*cos(omega*t);
nz = omega > 0;
if any(nz)
    alpha(nz, :) = alpha(nz, :) + (b(nz)./omega(nz)).*sin(omega(nz)*t);
end
if any(~nz)
    alpha(~nz, :) = a(~nz)*ones(1, numel(t)) + b(~nz)*t;
end
end

function out = alpha_of(a, b, omega, t, which)
t = t(:)';
if strcmp(which, 'val')
    out = modal_analytic(a, b, omega, t);
else
    out = zeros(numel(a), numel(t));
    nz  = omega > 0;
    out(nz, :)  = -(a(nz).*omega(nz)).*sin(omega(nz)*t) + b(nz).*cos(omega(nz)*t);
    out(~nz, :) = b(~nz)*ones(1, numel(t));
end
end

% ----------------------------------------------------------------------
function A = display_scale(W)
if size(W, 3) > 1, V = abs(W(:, :, 2:end)); else, V = abs(W); end
V = sort(V(:));
if isempty(V) || V(end) == 0, A = 1; return, end
A = V(max(1, ceil(0.995*numel(V))));
A = min(max(A, 0.02*V(end)), V(end));
end

% ----------------------------------------------------------------------
function opt = parse_opts(opt, args)
if numel(args) == 1 && isstruct(args{1})
    fn = fieldnames(args{1});
    for i = 1:numel(fn), opt = assign(opt, fn{i}, args{1}.(fn{i})); end
    return
end
if mod(numel(args), 2) ~= 0
    error('wave_solver_polygon:badArgs', 'Options must be Name/Value pairs.');
end
for i = 1:2:numel(args), opt = assign(opt, args{i}, args{i+1}); end
end

function opt = assign(opt, name, value)
fn = fieldnames(opt);
ix = find(strcmpi(fn, name), 1);
if isempty(ix)
    error('wave_solver_polygon:unknownOption', 'Unknown option ''%s''.', name);
end
opt.(fn{ix}) = value;
end
