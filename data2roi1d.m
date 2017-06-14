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
## @anchor{data2roi1d}
## @deftypefn {Function File} {@var{X} =} data2roi1d (@var{Xroi}, @var{ab}, @var{roi})
## Uses one-dimensional spectral data inside regions of interest from a data
## matrix to reconstruct portions of a full spectral dataset.
## @end deftypefn

function X = data2roi1d (Xroi, ab, roi)
  % check the number of input arguments.
  if (nargin != 3 || nargout != 1 || !ismatrix(Xroi) || !isvector(ab))
    % unknown. throw an exception.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 2 || rows(roi) < 1)
    % invalid. throw an exception.
    error('data2roi1d: invalid region of interest matrix');
  end

  % get the size of the output matrix.
  K = length(ab);
  N = rows(Xroi);

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % build the boundary matrix from the roi matrix.
  v = sort(findnearest(ab, roi), 2);

  % initialize the output matrix.
  X = zeros(N, K);

  % loop through the roi boundaries.
  k = 1;
  for i = 1 : rows(v)
    % get the roi start and end indices.
    k1 = v(i, 1);
    k2 = v(i, 2);

    % get the roi variable count.
    nk = k2 - k1 + 1;

    % extract the segment of interest.
    X(:, k1 : k2) = Xroi(:, k : k + nk - 1);

    % shift past the extracted region.
    k += nk;
  end
end

