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
## @anchor{binmanual2d}
## @deftypefn {Function File} {@var{xnew} =} binmanual2d (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binmanual2d (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binmanual2d (@var{X}, @var{ab}, @var{roi})
## Manually bin a two-dimensional spectrum or spectral dataset in @var{X} based
## on regions of interest provided in @var{roi}, or centers and widths provided
## in @var{centers} and @var{widths}. If regions of interest are used to bin,
## the @var{abnew} and @var{widths} values are optionally returnable.
## @end deftypefn

function [xnew, abnew, widths] = binmanual2d (X, ab, roi)
  % check for proper arguments.
  if (nargin != 3 || nargout != 3 || !iscell(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 4 || rows(roi) < 1)
    % invalid. throw an exception.
    error('binmanual2d: invalid region of interest matrix');
  end

  % get the size of the input matrix.
  [Ka, Kb] = size(X{1});
  N = length(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab{1}) != Kb || length(ab{2}) != Ka)
    % no match. throw an exception.
    error('binmanual2d: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a cell array of column vectors.
  ab{1} = reshape(ab{1}, Kb, 1);
  ab{2} = reshape(ab{2}, Ka, 1);

  % build the boundary matrix from the roi matrix.
  v1 = sort(findnearest(ab{1}, roi(:, [1, 2])), 2);
  v2 = sort(findnearest(ab{2}, roi(:, [3, 4])), 2);
  v = [v1, v2];

  % initialize output variables.
  xnew = [];
  abnew = [];
  widths = [];

  % loop through the bin boundaries.
  for idx = 1 : rows(v)
    % get the first dimension bin start and end indices.
    i1 = v(idx, 1);
    i2 = v(idx, 2);

    % get the second dimension bin start and end indices.
    j1 = v(idx, 3);
    j2 = v(idx, 4);

    % build the abscissas of the current bin.
    ab1 = ab{1}(i1 : i2);
    ab2 = ab{2}(j1 : j2);

    % build the next abscissa center, width and index.
    abnew = [abnew; median(ab1), median(ab2)];
    widths = [widths; abs(range(ab1)), abs(range(ab2))];

    % loop for each element of the cell array.
    for k = 1 : N
      % get the bin segment from the current matrix.
      Xk = X{k};
      Xi = Xk(j1 : j2, i1 : i2);

      % bin the segment into its row.
      xnew(k,idx) = abs(trapz(ab2, trapz(ab1, Xi, 2)));
    end
  end
end

