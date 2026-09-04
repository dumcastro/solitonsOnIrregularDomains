function test_wave_solver_polygon()
%TEST_WAVE_SOLVER_POLYGON  Verification of the Phase 2 conformal framework.
%
%   The conformal machinery has two exact cases to lean on.  Under the
%   identity map the weighted eigenproblem must collapse to the ordinary
%   Bessel spectrum of the unit disc, and under z = s*zeta it must give
%   that spectrum divided by s, because the image is a disc of radius s.
%   Both are checked to near machine precision.  For a genuinely non-radial
%   image the exact answer is unknown, so the tests fall back on internal
%   consistency: convergence of the lowest eigenvalue under refinement,
%   conservation of energy, and agreement between an analytic map
%   derivative and the finite-difference fallback.

fprintf('\n=== wave_solver_polygon verification ===\n\n');
R = struct('name', {}, 'err', {}, 'tol', {});
common = {'LambdaMax', 12, 'nT', 5, 'T', 1, 'Verbose', false};

% -- 1. identity map reproduces the disc spectrum -----------------------
s1 = wave_solver_polygon('MapFun', @(z) z, 'MapDiff', @(z) ones(size(z)), common{:});
E  = disc_eigenmodes(1, 'dirichlet', 12);
n  = min(20, numel(s1.mu));
R  = add(R, 'identity map -> disc spectrum', ...
         max(abs(s1.mu(1:n) - E.lambda(1:n))./E.lambda(1:n)), 1e-10);

% -- 2. a dilation scales the spectrum by 1/s ---------------------------
sc = 1.7;
s2 = wave_solver_polygon('MapFun', @(z) sc*z, 'MapDiff', @(z) sc*ones(size(z)), ...
                         common{:});
R  = add(R, 'dilation -> spectrum/s', ...
         max(abs(s2.mu(1:n) - E.lambda(1:n)/sc)./(E.lambda(1:n)/sc)), 1e-10);

% -- 3. the finite-difference fallback matches an analytic derivative ---
f  = @(z) z + 0.22*z.^3;
df = @(z) 1 + 0.66*z.^2;
sa = wave_solver_polygon('MapFun', f, 'MapDiff', df, common{:});
sn = wave_solver_polygon('MapFun', f, common{:});
R  = add(R, 'FD map derivative matches exact', ...
         max(abs(sa.mu(1:n) - sn.mu(1:n))./sa.mu(1:n)), 1e-6);

% -- 4. the lowest eigenvalues converge under refinement ----------------
lamList = [8 12 16 20];
mu1 = zeros(size(lamList)); mu4 = mu1;
for i = 1:numel(lamList)
    si = wave_solver_polygon('MapFun', f, 'MapDiff', df, 'LambdaMax', lamList(i), ...
                             'nT', 2, 'T', 1, 'Verbose', false);
    mu1(i) = si.mu(1); mu4(i) = si.mu(4);
end
R = add(R, 'mu_1 converges (successive diff)', abs(mu1(end) - mu1(end-1))/mu1(end), 1e-4);
R = add(R, 'mu_4 converges (successive diff)', abs(mu4(end) - mu4(end-1))/mu4(end), 1e-3);
R = add(R, 'mu_1 decreasing under refinement', max(diff(mu1)), 1e-12);

% -- 5. energy conservation on the mapped domain ------------------------
se = wave_solver_polygon('MapFun', f, 'MapDiff', df, 'LambdaMax', 14, ...
                         'T', 4, 'nT', 41, 'Verbose', false);
ee = se.energy(se.t);
R  = add(R, 'energy conserved on the polygon', max(abs(ee - ee(1)))/ee(1), 1e-12);

% -- 6. the conformal Jacobian really is the area element ---------------
% INT_D rho dA must equal the area enclosed by f(dD).  The reference comes
% from Green's theorem on the boundary alone, A = (1/2) INT Im(conj(z) z')
% with z(th) = f(exp(i th)) and z' obtained by Fourier differentiation, so
% it uses only values of f - never f' - and is spectrally accurate.
N  = 4096;
th = 2*pi*(0:N-1)/N;
zb = f(exp(1i*th));
kk = [0:N/2-1, -N/2:-1];
zp = ifft(1i*kk.*fft(zb));
area_boundary = 0.5*abs(sum(imag(conj(zb).*zp))*(2*pi/N));
[rq, wq] = gauss_legendre(400, 0, 1);
Nth = 512; tq = 2*pi*(0:Nth-1)/Nth;
ZT  = rq*exp(1i*tq);
area_quad = sum((wq.*rq)'*abs(df(ZT)).^2)*(2*pi/Nth);
R = add(R, 'conformal Jacobian gives the area', ...
        abs(area_quad - area_boundary)/area_boundary, 1e-12);

% -- report --------------------------------------------------------------
nfail = 0;
fprintf('  %-36s %12s %10s   %s\n', 'check', 'error', 'tol', 'result');
fprintf('  %s\n', repmat('-', 1, 74));
for i = 1:numel(R)
    ok = R(i).err <= R(i).tol;
    nfail = nfail + ~ok;
    if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-36s %12.3e %10.1e   %s\n', R(i).name, R(i).err, R(i).tol, verdict);
end
fprintf('  %s\n', repmat('-', 1, 74));
if nfail == 0
    fprintf('  all %d checks passed\n\n', numel(R));
else
    fprintf('  %d of %d checks FAILED\n\n', nfail, numel(R));
    error('test_wave_solver_polygon:failed', '%d check(s) failed.', nfail);
end
end

function R = add(R, name, err, tol)
R(end+1) = struct('name', name, 'err', err, 'tol', tol); %#ok<AGROW>
end
