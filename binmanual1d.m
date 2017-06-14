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
## @anchor{binmanual1d}
## @deftypefn {Function File} {@var{xnew} =} binmanual1d (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binmanual1d (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binmanual1d (@var{X}, @var{ab}, @var{roi})
## Manually bin a one-dimensional spectrum or spectral dataset in @var{X} based
## on regions of interest provided in @var{roi}, or centers and widths provided
## in @var{centers} and @var{widths}. If regions of interest are used to bin,
## the @var{abnew} and @var{widths} values are optionally returnable.
## @end deftypefn

function [xnew, abnew, widths] = binmanual1d (X, ab, roi)
  % check for proper arguments.
  if (nargin != 3 || nargout != 3 || !ismatrix(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 2 || rows(roi) < 1)
    % invalid. throw an exception.
    error('binmanual1d: invalid region of interest matrix');
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('binmanual1d: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % build the boundary matrix from the roi matrix.
  v = sort(findnearest(ab, roi), 2);

  % initialize output variables.
  xnew = [];
  abnew = [];
  widths = [];

  % loop through the bin boundaries.
  for i = 1 : rows(v)
    % get the bin start and end indices.
    k1 = v(i, 1);
    k2 = v(i, 2) - 1;

    % build the segments over which to calculate the next bin.
    Xi = X(:, k1 : k2);
    abi = ab(k1 : k2);

    % build the next output bin, abscissa center, width and index.
    xnew = [xnew, abs(trapz(abi, Xi, 2))];
    abnew = [abnew; median(abi)];
    widths = [widths; abs(range(abi))];
  end
end

