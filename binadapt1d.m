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
## @anchor{binadapt1d}
## @deftypefn {Function File} {@var{xnew} =} binadapt1d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binadapt1d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binadapt1d (@var{X}, @var{ab}, @var{w})
## Adaptively bin a one-dimensional spectrum or spectral dataset.
##
## It is highly recommended that you use @ref{binadapt} instead of this
## function directly.
## @end deftypefn

function [xnew, abnew, widths] = binadapt1d (X, ab, w, R)
  % check for proper arguments.
  if (!any(nargin == [2 : 4]) || nargout != 3 || !ismatrix(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the input matrix is actually a vector (single observation case).
  if (isvector(X))
    % reshape into a row vector to obey the matrix rules.
    X = reshape(X, 1, length(X));
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('binadapt: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % check if a third argument was supplied.
  if (nargin < 3 || isempty(w))
    % no. use the default initial bin width of 0.025 abscissa units.
    w = 0.025;
  end

  % check if a fourth argument was supplied.
  if (nargin < 4 || isempty(R))
    % no. use the default resolution parameter of 0.5.
    R = 0.5;
  end

  % initialize output variables.
  v = [];
  xnew = [];
  abnew = [];
  widths = [];

  % find the jump locations in the abscissa.
  kstops = findjumps(ab);
  kstops = [1; kstops; K + 1];

  % loop through the stops array.
  for i = 1 : length(kstops) - 1
    % get the sticky indices of the data (the ones that bridge the jump).
    k1 = kstops(i);
    k2 = kstops(i + 1) - 1;

    % slice the data matrix based on the current boundaries.
    Xi = X(:, k1 : k2);
    abi = ab(k1 : k2);

    % calculate the minimum bin size in points.
    kmin = floor(w * length(abi) / range(abi));

    % run the adaptive intelligent binning subroutine.
    vi = __binadapt1d(Xi, kmin, R) + k1;

    % append the results to the total boundary matrix.
    v = [v; vi];
  end

  % sort the boundary matrix.
  [jnk, idx] = sort(v(:, 1));
  v = v(idx, :);

  % loop through the bin boundaries.
  for i = 1 : rows(v)
    % get the bin start and end indices.
    k1 = v(i, 1);
    k2 = v(i, 2);

    % build the segments over which to calculate the next bin.
    Xi = X(:, k1 : k2);
    abi = ab(k1 : k2)';

    % build the next output bin, abscissa center, width and index.
    xnew = [xnew, abs(trapz(abi, Xi, 2))];
    abnew = [abnew; median(abi)];
    widths = [widths; abs(range(abi))];
  end
end

