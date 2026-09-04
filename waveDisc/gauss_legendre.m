function [x, w] = gauss_legendre(n, a, b)
%GAUSS_LEGENDRE  Gauss-Legendre nodes and weights (Golub-Welsch).
%
%   [X,W] = GAUSS_LEGENDRE(N) returns the N-point Gauss-Legendre rule on
%   [-1,1], so that SUM(W.*F(X)) approximates INTEGRAL(F,-1,1) exactly for
%   polynomials of degree <= 2N-1.
%
%   [X,W] = GAUSS_LEGENDRE(N,A,B) maps the rule to the interval [A,B].
%
%   X and W are returned as column vectors.  No toolbox required: the nodes
%   are the eigenvalues of the symmetric tridiagonal Jacobi matrix and the
%   weights come from the first components of its eigenvectors.

if n < 1
    error('gauss_legendre:badN', 'N must be a positive integer.');
end

if n == 1
    x = 0; w = 2;
else
    k    = (1:n-1)';
    beta = k./sqrt(4*k.^2 - 1);          % off-diagonal of the Jacobi matrix
    T    = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    [x, ix] = sort(diag(D));
    w = 2*(V(1, ix).^2)';
end

x = x(:); w = w(:);

if nargin > 1
    if nargin < 3
        error('gauss_legendre:badInterval', 'Supply both A and B, or neither.');
    end
    x = 0.5*(b - a)*x + 0.5*(a + b);
    w = 0.5*(b - a)*w;
end
