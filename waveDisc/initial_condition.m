function [w0, v0, info] = initial_condition(name, opt)
%INITIAL_CONDITION  Smooth (C^2) initial data for the disc wave solver.
%
%   [W0,V0] = INITIAL_CONDITION(NAME,OPT) returns function handles
%   W0(X,Y) and V0(X,Y) giving the initial displacement and the initial
%   velocity, both accepting array arguments elementwise.
%
%   OPT is the option struct built by WAVE_SOLVER_DISC; the fields each
%   profile reads are listed below.  NAME is one of
%
%     'gaussian'   (default)  a single Gaussian bump, the standard option
%                  w0 = A exp(-((x-x0)^2 + (y-y0)^2)/(2 sigma^2))
%                  fields: amplitude, x0, y0, sigma
%
%     'ring'       a radially symmetric Gaussian annulus
%                  w0 = A exp(-(r - r0)^2/(2 sigma^2))
%                  fields: amplitude, r0, sigma
%
%     'twobumps'   two Gaussians of opposite sign, for interference
%                  fields: amplitude, x0, y0, sigma
%
%     'mode'       a single exact eigenfunction J_m(lambda r) cos/sin(m th),
%                  for which the solution is known in closed form and which
%                  is therefore the natural verification case
%                  fields: modeM, modeN, modePar, R, BC, amplitude
%
%     'custom'     user-supplied handles
%                  fields: w0fun (required), v0fun (optional)
%
%   Every profile takes an optional initial velocity through OPT.v0fun; if
%   that is empty the field starts from rest, V0 = 0, which is the classical
%   "plucked membrane" setup.
%
%   The third output INFO is a struct describing the profile: INFO.label is
%   a short human-readable name and INFO.exact, when present, is a handle
%   EXACT(X,Y,T) giving the closed-form solution (only for 'mode').
%
%   Note on smoothness: a Gaussian does not vanish identically on r = R, so
%   under Dirichlet conditions the eigenfunction expansion converges slowly
%   unless the bump is well inside the disc.  WAVE_SOLVER_DISC measures the
%   residual max|w0| on the boundary and warns when it is not negligible.
%
%   See also WAVE_SOLVER_DISC, DISC_EIGENMODES.

if nargin < 1 || isempty(name), name = 'gaussian'; end

info = struct('label', name);

switch lower(name)

    case 'gaussian'
        A = opt.amplitude; x0 = opt.x0; y0 = opt.y0; s = opt.sigma;
        w0 = @(x, y) A*exp(-((x - x0).^2 + (y - y0).^2)/(2*s^2));
        info.label = sprintf('Gaussian, centre (%.2f, %.2f), \\sigma = %.3g', x0, y0, s);

    case 'ring'
        A = opt.amplitude; r0 = opt.r0; s = opt.sigma;
        w0 = @(x, y) A*exp(-(sqrt(x.^2 + y.^2) - r0).^2/(2*s^2));
        info.label = sprintf('Gaussian ring, r_0 = %.3g, \\sigma = %.3g', r0, s);

    case 'twobumps'
        A = opt.amplitude; x0 = opt.x0; y0 = opt.y0; s = opt.sigma;
        w0 = @(x, y) A*exp(-((x - x0).^2 + (y - y0).^2)/(2*s^2)) ...
                   - A*exp(-((x + x0).^2 + (y + y0).^2)/(2*s^2));
        info.label = sprintf('Two opposed Gaussians at (\\pm%.2f, \\pm%.2f), \\sigma = %.3g', x0, y0, s);

    case 'mode'
        m = opt.modeM; n = opt.modeN; A = opt.amplitude;
        if strncmpi(opt.BC, 'n', 1), kind = 'Jprime'; else, kind = 'J'; end
        z   = bessel_zeros(m, n, kind);
        j   = z(n);
        lam = j/opt.R;
        if opt.modePar >= 0
            ang = @(th) cos(m*th);  ps = 'cos';
        else
            ang = @(th) sin(m*th);  ps = 'sin';
        end
        w0 = @(x, y) A*besselj(m, lam*sqrt(x.^2 + y.^2)).*ang(atan2(y, x));
        info.label  = sprintf('eigenmode (m,n) = (%d,%d), %s, \\lambda = %.4f', m, n, ps, lam);
        info.lambda = lam;
        info.exact  = @(x, y, t) w0(x, y).*cos(opt.c*lam*t);   % v0 = 0

    case 'custom'
        if ~isfield(opt, 'w0fun') || isempty(opt.w0fun)
            error('initial_condition:noHandle', ...
                'IC ''custom'' requires the option ''w0fun'', a handle @(x,y).');
        end
        w0 = opt.w0fun;
        info.label = 'custom w_0(x,y)';

    otherwise
        error('initial_condition:badName', ...
            ['Unknown initial condition ''%s''. Choose gaussian | ring | ' ...
             'twobumps | mode | custom.'], name);
end

if isfield(opt, 'v0fun') && ~isempty(opt.v0fun)
    v0 = opt.v0fun;
    info.label = [info.label, ', with initial velocity'];
else
    v0 = @(x, y) zeros(size(x));
end
