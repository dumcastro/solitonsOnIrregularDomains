function test_wave_solver_disc()
%TEST_WAVE_SOLVER_DISC  Verification suite for the disc wave solver.
%
%   Runs a sequence of independent checks - exact Bessel zeros, eigenmode
%   normalisation, reproduction of the initial data, closed-form single-mode
%   evolution, boundary conditions, energy conservation, agreement between
%   the analytic and ODE45 time integrators, conservation of mass under
%   Neumann conditions, preservation of radial symmetry, and a direct
%   finite-difference residual of the PDE itself - and prints a PASS/FAIL
%   table.  Errors are raised only at the end, so a full report is always
%   produced.
%
%   Usage:  test_wave_solver_disc

fprintf('\n=== wave_solver_disc verification ===\n\n');
R = struct('name', {}, 'err', {}, 'tol', {});

% -- 1. Bessel zeros against tabulated values -------------------------
ref = [2.404825557695773, 5.520078110286311, 8.653727912911013];
R = add(R, 'Bessel zeros j_{0,n}', max(abs(bessel_zeros(0,3) - ref(:))), 1e-12);
refp = [1.841183781340659, 5.331442773525032, 8.536316366346286];
R = add(R, 'Bessel zeros j''_{1,n}', max(abs(bessel_zeros(1,3,'Jprime') - refp(:))), 1e-12);

% -- 2. eigenmode norms by brute-force quadrature ----------------------
E = disc_eigenmodes(1, 'dirichlet', 12);
[rr, wr] = gauss_legendre(400, 0, 1);
Nth = 512; th = 2*pi*(0:Nth-1)/Nth;
e = 0;
for k = 1:E.K
    if E.par(k) > 0, ang = cos(E.m(k)*th); else, ang = sin(E.m(k)*th); end
    P = besselj(E.m(k), E.lambda(k)*rr)*ang;
    num = sum((wr.*rr)'*(P.^2))*(2*pi/Nth);
    e = max(e, abs(num - E.nrm(k))/E.nrm(k));
end
R = add(R, 'eigenmode norms (Dirichlet)', e, 1e-10);

EN = disc_eigenmodes(1, 'neumann', 12);
e = 0;
for k = 1:EN.K
    if EN.par(k) > 0, ang = cos(EN.m(k)*th); else, ang = sin(EN.m(k)*th); end
    P = besselj(EN.m(k), EN.lambda(k)*rr)*ang;
    num = sum((wr.*rr)'*(P.^2))*(2*pi/Nth);
    e = max(e, abs(num - EN.nrm(k))/EN.nrm(k));
end
R = add(R, 'eigenmode norms (Neumann)', e, 1e-10);

% -- 3. reproduction of the initial data at t = 0 ----------------------
s = wave_solver_disc('sigma', 0.12, 'x0', 0.2, 'T', 3, 'nT', 11, ...
                     'Verbose', false, 'PrecomputeField', false);
[xg, yg] = meshgrid(linspace(-0.95, 0.95, 25));
in = hypot(xg, yg) <= 0.95;
w0n = s.eval(xg(in), yg(in), 0);
w0e = s.w0fun(xg(in), yg(in));
R = add(R, 'w(.,0) reproduces w0', max(abs(w0n - w0e)), 1e-8);

% -- 4. single eigenmode against its closed form -----------------------
sm = wave_solver_disc('IC', 'mode', 'modeM', 3, 'modeN', 4, 'c', 1.3, ...
                      'T', 5, 'nT', 21, 'Verbose', false, 'PrecomputeField', false);
tt = sm.t;
Wn = sm.eval(xg(in), yg(in), tt);
We = sm.icInfo.exact(xg(in), yg(in), tt);
R = add(R, 'single mode vs closed form', max(abs(Wn(:) - We(:))), 1e-10);

% -- 5. boundary conditions --------------------------------------------
thb = linspace(0, 2*pi, 120)';
Wb = s.eval(cos(thb), sin(thb), s.t);
R = add(R, 'Dirichlet w = 0 on dD', max(abs(Wb(:))), 1e-12);

sn = wave_solver_disc('BC', 'neumann', 'sigma', 0.12, 'x0', 0.2, 'T', 3, ...
                      'nT', 11, 'Verbose', false, 'PrecomputeField', false);
h  = 1e-4;
dW = (sn.eval((1-2*h)*cos(thb), (1-2*h)*sin(thb), sn.t) ...
      - 4*sn.eval((1-h)*cos(thb), (1-h)*sin(thb), sn.t) ...
      + 3*sn.eval(cos(thb), sin(thb), sn.t))/(2*h);
R = add(R, 'Neumann dw/dn = 0 on dD', max(abs(dW(:))), 1e-6);

% -- 6. energy conservation ---------------------------------------------
Ee = s.energy(s.t);
R = add(R, 'energy conserved (analytic)', max(abs(Ee - Ee(1)))/Ee(1), 1e-13);

% -- 7. mass conservation under Neumann ---------------------------------
% With dw/dn = 0 and w_t(.,0) = 0 the mean of w over the disc is constant.
[rr2, wr2] = gauss_legendre(60, 0, 1);
Nth2 = 64; th2 = 2*pi*(0:Nth2-1)/Nth2;
Xm = rr2*cos(th2); Ym = rr2*sin(th2);
mass = zeros(1, numel(sn.t));
for i = 1:numel(sn.t)
    Wi = sn.eval(Xm, Ym, sn.t(i));
    mass(i) = sum((wr2.*rr2)'*Wi)*(2*pi/Nth2);
end
R = add(R, 'Neumann mass conserved', max(abs(mass - mass(1)))/abs(mass(1)), 1e-12);

% -- 8. radial symmetry of a centred pulse ------------------------------
sc = wave_solver_disc('x0', 0, 'y0', 0, 'sigma', 0.2, 'T', 2, 'nT', 11, ...
                      'Verbose', false, 'PrecomputeField', false);
rtest = 0.37;
Wr = sc.eval(rtest*cos(thb), rtest*sin(thb), sc.t);
R = add(R, 'centred pulse stays radial', max(max(Wr, [], 1) - min(Wr, [], 1)), 1e-12);

% -- 9. ODE45 agrees with the closed form -------------------------------
ode45args = {'IC', 'mode', 'modeM', 2, 'modeN', 3, 'T', 4, 'nT', 9, ...
             'Verbose', false, 'PrecomputeField', false};
so = wave_solver_disc(ode45args{:}, 'TimeScheme', 'ode45');
sa = wave_solver_disc(ode45args{:});
R = add(R, 'ode45 vs analytic (modal)', ...
        max(abs(so.alpha(:) - sa.alpha(:)))/max(abs(sa.alpha(:))), 1e-7);

% -- 10. direct finite-difference residual of the PDE -------------------
% Checks w_tt - c^2 (w_xx + w_yy) = 0 at interior points, using the solver
% only as a black box.  This is independent of how the solver is built.
sp = wave_solver_disc('sigma', 0.14, 'x0', 0.2, 'c', 1.0, 'T', 2, 'nT', 5, ...
                      'Verbose', false, 'PrecomputeField', false);
hh = 2e-3; kk = 2e-3; t0 = 0.8;
pts = [0.10 0.15; -0.30 0.20; 0.05 -0.40; 0.45 0.10];
res = 0; scale = 0;
for p = 1:size(pts, 1)
    x = pts(p, 1); y = pts(p, 2);
    wtt = (sp.eval(x, y, t0+kk) - 2*sp.eval(x, y, t0) + sp.eval(x, y, t0-kk))/kk^2;
    wxx = (sp.eval(x+hh, y, t0) - 2*sp.eval(x, y, t0) + sp.eval(x-hh, y, t0))/hh^2;
    wyy = (sp.eval(x, y+hh, t0) - 2*sp.eval(x, y, t0) + sp.eval(x, y-hh, t0))/hh^2;
    res   = max(res, abs(wtt - sp.opt.c^2*(wxx + wyy)));
    scale = max(scale, abs(wtt));
end
R = add(R, 'PDE residual w_tt - c^2 Lap w', res/scale, 5e-5);

% -- 11. spectral convergence of the projection -------------------------
% For data that is smooth and safely interior, widening the cutoff must
% drive the projection error down super-algebraically, until it saturates
% at the level by which the Gaussian tail violates the boundary condition.
lamList = [12 18 24 30 36];
pe = zeros(size(lamList));
for i = 1:numel(lamList)
    si = wave_solver_disc('sigma', 0.15, 'x0', 0.15, 'LambdaMax', lamList(i), ...
                          'nT', 2, 'Verbose', false, 'PrecomputeField', false);
    pe(i) = si.diag.projErr;
end
R = add(R, 'projection error monotone', max(diff(pe)), 0);
R = add(R, 'projection error spectral', pe(end)/pe(1), 1e-5);

% -- report --------------------------------------------------------------
nfail = 0;
fprintf('  %-34s %12s %10s   %s\n', 'check', 'error', 'tol', 'result');
fprintf('  %s\n', repmat('-', 1, 72));
for i = 1:numel(R)
    ok = R(i).err <= R(i).tol;
    nfail = nfail + ~ok;
    if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-34s %12.3e %10.1e   %s\n', R(i).name, R(i).err, R(i).tol, verdict);
end
fprintf('  %s\n', repmat('-', 1, 72));
if nfail == 0
    fprintf('  all %d checks passed\n\n', numel(R));
else
    fprintf('  %d of %d checks FAILED\n\n', nfail, numel(R));
    error('test_wave_solver_disc:failed', '%d check(s) failed.', nfail);
end
end

function R = add(R, name, err, tol)
R(end+1) = struct('name', name, 'err', err, 'tol', tol); %#ok<AGROW>
end
