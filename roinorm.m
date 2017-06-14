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
## @anchor{roinorm}
## @deftypefn {Function File} {@var{Xn} =} roinorm (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{Xn}, @var{s}] =} roinorm (@var{X}, @var{ab}, @var{roi})
## Normalize the observations of a data matrix to such that the integral
## of a spectral region specified by @var{roi} is one. The calculated
## normalization factors may be optionally returned in @var{s}. The
## values specified in @var{roi} must correspond to those in the abscissa
## @var{ab}.
## @end deftypefn

function [Xnew, s] = roinorm (X, ab, roi)
  % check for proper arguments.
  if (nargin != 3 || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the input matrix is complex.
  if (iscomplex(X))
    % only use the real portion.
    X = real(X);
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check the abscissa vector.
  if (!isvector(ab) || length(ab) != K)
    % invalid. throw an exception.
    error('roinorm: variable count and abscissa length do not match');
  end

  % check the roi vector.
  if (!isvector(roi) || length(roi) != 2)
    % invalid. throw an exception.
    error('roinorm: invalid region of interest vector');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % reshape the roi vector into a row vector.
  roi = reshape(roi, 1, 2);

  % build the boundary vector from the roi matrix.
  v = sort(findnearest(ab, roi), 2);

  % calculate the integrals of the regions of interest.
  ints = integrals(X, ab, roi);

  % divide the rows of the data matrix by its own roi integrals.
  Xnew = diag(1 ./ ints) * X;

  % check if a second argument was requested.
  if (nargout == 2)
    % yes. return the normalization factors.
    s = ints;
  end
end

