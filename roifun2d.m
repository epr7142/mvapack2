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
## @anchor{roifun2d}
## @deftypefn {Function File} {@var{Y} =} roifun2d (@var{X}, @var{ab}, @var{roi}, @var{func})
## Performs a function on two-dimensional spectral data inside regions
## of interest.
## @end deftypefn

function Y = roifun2d (X, ab, roi, func)
  % check the number of input arguments.
  if (nargin != 4 || nargout != 1 || !iscell(X) || !isvector(ab))
    % unknown. throw an exception.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 4 || rows(roi) < 1)
    % invalid. throw an exception.
    error('roifun2d: invalid region of interest matrix');
  end

  % get the size of the input matrix.
  [Ka, Kb] = size(X{1});
  N = length(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab{1}) != Kb || length(ab{2}) != Ka)
    % no match. throw an exception.
    error('roifun2d: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab{1} = reshape(ab{1}, Kb, 1);
  ab{2} = reshape(ab{2}, Ka, 1);

  % build the boundary matrix from the roi matrix.
  v1 = sort(findnearest(ab{1}, roi(:, [1, 2])), 2);
  v2 = sort(findnearest(ab{2}, roi(:, [3, 4])), 2);
  v = [v1, v2];

  % initialize the output variable.
  Y = [];

  % loop through the roi boundaries.
  for idx = 1 : rows(v)
    % get the first dimension bin start and end indices.
    i1 = v(idx, 1);
    i2 = v(idx, 2);

    % get the second dimension bin start and end indices.
    j1 = v(idx, 3);
    j2 = v(idx, 4);

    % get the bin sizes.
    ni = i2 - i1 + 1;
    nj = j2 - j1 + 1;
    nk = ni * nj;

    % build the abscissas of the current roi.
    ab1 = ab{1}(i1 : i2);
    ab2 = ab{2}(j1 : j2);

    % initialize the data matrix segment.
    Xi = cell(N, 1);

    % build the next abscissa segment.
    abi = [median(ab1), median(ab2)];

    % loop for each element of the cell array.
    for k = 1 : N
      % get the roi segment from the current matrix.
      Xk = X{k};
      Xi{k} = Xk(j1 : j2, i1 : i2);
    end

    % append the extracted segment into the outputs.
    Y = [Y, func(Xi, abi)];
  end
end

