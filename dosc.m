## Copyright (C) 2013 University of Nebraska Board of Regents.
## Written by Bradley Worley <bradley.worley@huskers.unl.edu>.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @anchor{dosc}
## @deftypefn {Function File} {[@var{Z}, @var{W}, @var{P}, @var{T}] =} dosc (@var{X}, @var{Y}, @var{A})
## @deftypefnx {Function File} {[@var{Z}, @var{W}, @var{P}, @var{T}] =} dosc (@var{X}, @var{Y}, @var{A}, @var{tol})
## Performs Direct Orthogonal Signal Correction as described in:
##
## @quotation
## Westerhuis J. A., de Jong S., Smilde A. K., `Direct orthogonal signal
## correction', Chemometrics and Intelligent Laboratory Systems,
## 56 (2001): 13-25.
## @end quotation
##
## The function requires a data matrix @var{X}, a response matrix @var{Y},
## and a number of OSC components to calculate @var{A}. If a fourth argument
## @var{tol} is provided, it will be used as the tolerance for calculating
## a pseudoinverse of the data matrix. Otherwise, a default value of 1e-3
## will be used.
##
## The function returns a corrected data matrix @var{Z}, OSC weights @var{W},
## OSC loadings @var{P} and OSC score components @var{T}.
## @end deftypefn

function [Z, W, P, T] = dosc (X, Y, A, tol)
  % check that the tolerance was provided.
  if (nargin < 4 || isempty(tol))
    % not provided. use the default value.
    tol = 1.0e-3;
  end

  Yhat = X * (pinv(X')' * Y);
  Xhat = X - Yhat * (pinv(Yhat) * X);

  [Txh, D] = eigs(Xhat * Xhat', A);
  pinvX = pinv(X', tol)';

  W = pinvX * Txh;
  T = X * W;
  P = X' * T * inv(T' * T);
  Z = X - T * P';
end

