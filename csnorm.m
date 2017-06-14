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
## @anchor{csnorm}
## @deftypefn {Function File} {@var{Xn} =} csnorm (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{s}] =} csnorm (@var{X})
## Normalize the observations of a data matrix to constant integral. The
## calculated normalization factors may be optionally returned in @var{s}.
## @end deftypefn

function [Xnew, s] = csnorm (X)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the input matrix is complex.
  if (iscomplex(X))
    % only use the real portion.
    X = real(X);
  end

  % divide the rows of the data matrix by its own row integrals.
  Xnew = diag(1 ./ trapz(X, 2)) * X;

  % check if a second argument was requested.
  if (nargout == 2)
    % yes. return the normalization factors.
    s = trapz(X, 2);
  end
end

