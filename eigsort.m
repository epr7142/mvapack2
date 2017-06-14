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
## @anchor{eigsort}
## @deftypefn {Function File} {[@var{V}, @var{lambda}] =} eigsort (@var{A})
## @deftypefnx {Function File} {[@var{V}, @var{lambda}] =} eigsort (@var{A}, @var{B})
## Calculates the eigendecomposition of @var{A} such that the eigenvectors and
## eigenvalues are sorted according to decreasing eigenvalue magnitude.
## Alternatively, a second argument may be supplied such that the result
## is an eigendecomposition of @math{B^{-1} A}.
## @end deftypefn

function [V, lambda] = eigsort (A, B)
  % check the arguments.
  if (nargin == 2)
    % run an eigendecomposition of inv(B)*A into V and L.
    [V, L] = eig(A, B);
  else
    % run a complete eigendecomposition of A into V and L.
    [V, L] = eig(A);
  end

  % get the diagonal elements of the eigenvalue matrix.
  lambda = diag(L);

  % sort the eigenvalues according to absolute value (e.g. if complex).
  [jnk, idx] = sort(abs(lambda));

  % rearrange the eigenvector matrix and eigenvalue vector such that the
  % entries follow the sorting order determined in the above sort step.
  V = V(:,flipud(idx));
  lambda = lambda(flipud(idx));
end

