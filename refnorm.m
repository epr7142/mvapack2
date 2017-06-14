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
## @anchor{refnorm}
## @deftypefn {Function File} {@var{Xn} =} refnorm (@var{X}, @var{ab}, @var{refcs})
## @deftypefnx {Function File} {[@var{Xn}, @var{s}] =} refnorm (@var{X}, @var{ab}, @var{refcs})
## Normalize the observations of a data matrix to such that the maximum
## intensity of a spectral region centered around @var{refcs} is one. The
## calculated normalization factors may be optionally returned in @var{s}.
## The values specified in @var{roi} must correspond to those in the abscissa
## @var{ab}.
## @end deftypefn

function [Xnew, s] = refnorm (X, ab, refcs, fuzz)
  % check for proper arguments.
  if (!any(nargin == [2 : 4]) || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a reference value was provided.
  if (nargin < 3 || isempty(refcs))
    % use the default value.
    refcs = 0.0;
  end

  % check if a fuzz value was provided.
  if (nargin < 4 || isempty(fuzz))
    % use the default value.
    fuzz = 10;
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
    error('refnorm: variable count and abscissa length do not match');
  end

  % check the reference value.
  if (!isscalar(refcs) || refcs < min(ab) || refcs > max(ab))
    % invalid. throw an exception.
    error('refnorm: invalid reference value');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % find the abscissa point index nearest to the reference value.
  refidx = findnearest(ab, refcs);

  % compute the normalization factors.
  ints = max(X(:, refidx - fuzz : refidx + fuzz), [], 2);

  % divide the rows of the data matrix by its own roi integrals.
  Xnew = diag(1 ./ ints) * X;

  % check if a second argument was requested.
  if (nargout == 2)
    % yes. return the normalization factors.
    s = ints;
  end
end

