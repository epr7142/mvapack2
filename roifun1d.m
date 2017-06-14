## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{roifun1d}
## @deftypefn {Function File} {@var{Y} =} roifun1d (@var{X}, @var{ab}, @var{roi}, @var{func})
## Performs a function on one-dimensional spectral data inside regions of
## interest.
## @end deftypefn

function Y = roifun1d (X, ab, roi, func)
  % check the number of input arguments.
  if (nargin != 4 || nargout != 1 || !ismatrix(X) || !isvector(ab))
    % unknown. throw an exception.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 2 || rows(roi) < 1)
    % invalid. throw an exception.
    error('roifun1d: invalid region of interest matrix');
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('roifun1d: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % build the boundary matrix from the roi matrix.
  v = sort(findnearest(ab, roi), 2);

  % initialize the output variable.
  Y = [];

  % loop through the roi boundaries.
  for i = 1 : rows(v)
    % get the roi start and end indices.
    k1 = v(i, 1);
    k2 = v(i, 2);

    % extract the segment of interest.
    Xi = X(:, k1 : k2);
    abi = ab(k1 : k2);

    % append the extracted segment to the outputs.
    Y = [Y, func(Xi, abi)];
  end
end

