function E = disc_eigenmodes(R, bc, lambdaMax, maxModes)
%DISC_EIGENMODES  Dirichlet/Neumann eigenbasis of -Laplacian on a disc.
%
%   E = DISC_EIGENMODES(R,BC,LAMBDAMAX) assembles every eigenfunction of
%
%       -Lap phi = lambda^2 phi   in  {r < R},   phi = 0      (Dirichlet)
%                                                dphi/dn = 0  (Neumann)
%
%   with lambda <= LAMBDAMAX.  In polar coordinates the eigenfunctions are
%
%       phi_{m,n}^cos = J_m(lambda_{mn} r) cos(m theta),   m >= 0
%       phi_{m,n}^sin = J_m(lambda_{mn} r) sin(m theta),   m >= 1
%
%   where lambda_{mn} = j_{m,n}/R with j_{m,n} the n-th positive zero of
%   J_m (Dirichlet) or of J_m' (Neumann).  The Neumann basis additionally
%   contains the constant mode lambda = 0.
%
%   BC is 'dirichlet' (default) or 'neumann'.  MAXMODES caps the basis size
%   (default 20000) and triggers an error if LAMBDAMAX asks for more.
%
%   The returned struct E has fields, all K-by-1 and sorted by lambda:
%       E.m       angular order
%       E.n       radial index
%       E.par     1 for the cos family, -1 for the sin family
%       E.lambda  eigenvalue sqrt, so that -Lap phi = lambda^2 phi
%       E.jz      the underlying Bessel (or Bessel-derivative) zero
%       E.nrm     squared L2(disc) norm  INT |phi|^2 r dr dtheta
%   plus E.R, E.bc, E.lambdaMax, E.K (basis size) and E.Mmax.
%
%   The squared norms are the closed forms
%       Dirichlet: nrm = C_m (R^2/2) J_{m+1}(j)^2
%       Neumann  : nrm = C_m (R^2/2) (1 - m^2/j^2) J_m(j)^2
%   with C_m = 2*pi for m = 0 and pi otherwise; the Neumann constant mode
%   has nrm = pi R^2.
%
%   See also BESSEL_ZEROS, WAVE_SOLVER_DISC.

if nargin < 2 || isempty(bc),        bc        = 'dirichlet'; end
if nargin < 3 || isempty(lambdaMax), lambdaMax = 40/R;        end
if nargin < 4 || isempty(maxModes),  maxModes  = 20000;       end

bc = lower(bc);
switch bc
    case {'dirichlet', 'd'}, bc = 'dirichlet'; kind = 'J';
    case {'neumann', 'n'},   bc = 'neumann';   kind = 'Jprime';
    otherwise
        error('disc_eigenmodes:badBC', 'BC must be ''dirichlet'' or ''neumann''.');
end

jMax = lambdaMax*R;         % cutoff expressed on the Bessel-zero axis

m_ = []; n_ = []; p_ = []; j_ = []; nr_ = [];

% Neumann: the constant eigenfunction phi = 1, lambda = 0.
if strcmp(bc, 'neumann')
    m_(end+1,1) = 0; n_(end+1,1) = 0; p_(end+1,1) = 1;
    j_(end+1,1) = 0; nr_(end+1,1) = pi*R^2;
end

m = 0;
while true
    % McMahon: j_{m,n} ~ (n + m/2 - 1/4)*pi, inverted to size the request.
    nz = max(1, ceil(jMax/pi - m/2 + 1/4) + 3);
    z  = bessel_zeros(m, nz, kind);
    z  = z(z <= jMax);

    if isempty(z)
        % j_{m,1} grows monotonically in m, so no larger m can contribute.
        break
    end

    Cm = 2*pi; if m >= 1, Cm = pi; end
    if strcmp(bc, 'dirichlet')
        rad = 0.5*R^2 * besselj(m+1, z).^2;
    else
        rad = 0.5*R^2 * (1 - m^2./z.^2) .* besselj(m, z).^2;
    end
    nrm = Cm*rad;

    nn = (1:numel(z))';
    m_  = [m_;  repmat(m, numel(z), 1)];   %#ok<AGROW>
    n_  = [n_;  nn];                       %#ok<AGROW>
    p_  = [p_;  ones(numel(z), 1)];        %#ok<AGROW>
    j_  = [j_;  z];                        %#ok<AGROW>
    nr_ = [nr_; nrm];                      %#ok<AGROW>

    if m >= 1                              % the sin partner, same lambda
        m_  = [m_;  repmat(m, numel(z), 1)];  %#ok<AGROW>
        n_  = [n_;  nn];                      %#ok<AGROW>
        p_  = [p_; -ones(numel(z), 1)];       %#ok<AGROW>
        j_  = [j_;  z];                       %#ok<AGROW>
        nr_ = [nr_; nrm];                     %#ok<AGROW>
    end

    if numel(m_) > maxModes
        error('disc_eigenmodes:tooManyModes', ...
            ['LambdaMax = %.3g needs more than MaxModes = %d basis functions.\n' ...
             'Lower LambdaMax (a wider initial pulse needs fewer modes) or ' ...
             'raise MaxModes.'], lambdaMax, maxModes);
    end
    m = m + 1;
end

lam = j_/R;
[~, ix] = sort(lam + 1e-12*double(p_ < 0));   % lambda order, cos before sin

E = struct();
E.R         = R;
E.bc        = bc;
E.lambdaMax = lambdaMax;
E.m      = m_(ix);
E.n      = n_(ix);
E.par    = p_(ix);
E.jz     = j_(ix);
E.lambda = lam(ix);
E.nrm    = nr_(ix);
E.K      = numel(ix);
E.Mmax   = max(E.m);
