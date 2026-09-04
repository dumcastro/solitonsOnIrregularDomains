function y = dbesselj(m, x)
%DBESSELJ  Derivative of the Bessel function of the first kind, dJ_m/dx.
%
%   Y = DBESSELJ(M,X) returns J_m'(X) using the standard recurrence
%       J_m'(x) = ( J_{m-1}(x) - J_{m+1}(x) ) / 2,      m >= 1
%       J_0'(x) = -J_1(x).
%
%   See also BESSELJ, BESSEL_ZEROS.

if m == 0
    y = -besselj(1, x);
else
    y = 0.5*(besselj(m-1, x) - besselj(m+1, x));
end
