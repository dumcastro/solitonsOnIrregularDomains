function z = bessel_zeros(m, nz, kind)
%BESSEL_ZEROS  First NZ positive zeros of J_m or of J_m'.
%
%   Z = BESSEL_ZEROS(M,NZ)            returns j_{M,1} .. j_{M,NZ},
%                                     the positive zeros of besselj(M,.).
%   Z = BESSEL_ZEROS(M,NZ,'Jprime')   returns j'_{M,1} .. j'_{M,NZ},
%                                     the positive zeros of J_M'.
%
%   The zeros are found by bracketing sign changes on a fine scan and
%   polishing each bracket with FZERO, which is robust for every order M
%   (unlike pure McMahon asymptotics, which degrade for small N).  Results
%   are cached per (KIND,M) between calls.
%
%   Note both J_M and J_M' are strictly positive on (0,M] for M >= 1, so the
%   scan may safely start at X = M and never mistakes the M-fold zero of
%   J_M at the origin for a simple zero.
%
%   See also BESSELJ, DBESSELJ, DISC_EIGENMODES.

if nargin < 3 || isempty(kind), kind = 'J'; end

switch lower(kind)
    case 'j',      f = @(x) besselj(m, x);
    case 'jprime', f = @(x) dbesselj(m, x);
    otherwise
        error('bessel_zeros:badKind', 'KIND must be ''J'' or ''Jprime''.');
end

% ---- cache ------------------------------------------------------------
persistent CACHE
if isempty(CACHE), CACHE = struct(); end
key = sprintf('%s_%d', lower(kind), m);
if isfield(CACHE, key)
    zc = CACHE.(key);
    if numel(zc) >= nz
        z = zc(1:nz);
        return
    end
end

% ---- scan for sign changes, polish with fzero -------------------------
dx = 0.05;                  % consecutive zeros are separated by ~pi
a  = max(0.05, m);          % no positive zero of J_m or J_m' lies in (0,m]
z  = zeros(0, 1);
chunk = max(40, (nz + 1)*pi + 10);

while numel(z) < nz
    b  = a + chunk;
    xs = (a:dx:b)';
    fs = f(xs);
    sgn = sign(fs);
    idx = find(sgn(1:end-1).*sgn(2:end) < 0);
    for k = 1:numel(idx)
        i = idx(k);
        z(end+1, 1) = fzero(f, [xs(i) xs(i+1)]); %#ok<AGROW>
    end
    hit = find(fs == 0);     % a scan node landing exactly on a zero
    if ~isempty(hit)
        z = [z; xs(hit)];    %#ok<AGROW>
    end
    a = b;
end

z = sort(z);
z = z([true; diff(z) > 1e-8]);   % drop duplicates from the two detections
z = z(1:nz);

CACHE.(key) = z;
