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
## @anchor{data2roi2d}
## @deftypefn {Function File} {@var{X} =} data2roi2d (@var{Xroi}, @var{ab}, @var{roi})
## Uses two-dimensional spectral data inside regions of interest from a data
## matrix to reconstruct portions of a full spectral dataset.
## @end deftypefn

function X = data2roi2d (Xroi, ab, roi)
  % check the number of input arguments.
  if (nargin != 3 || nargout != 1 || !ismatrix(Xroi) || !isvector(ab))
    % unknown. throw an exception.
    print_usage();
  end

  % check the roi matrix.
  if (!ismatrix(roi) || columns(roi) != 4 || rows(roi) < 1)
    % invalid. throw an exception.
    error('data2roi2d: invalid region of interest matrix');
  end

  % get the size of the output matrix.
  Ka = length(ab{2});
  Kb = length(ab{1});
  N = rows(Xroi);

  % reshape the abscissa into a column vector.
  ab{1} = reshape(ab{1}, Kb, 1);
  ab{2} = reshape(ab{2}, Ka, 1);

  % build the boundary matrix from the roi matrix.
  v1 = sort(findnearest(ab{1}, roi(:, [1, 2])), 2);
  v2 = sort(findnearest(ab{2}, roi(:, [3, 4])), 2);
  v = [v1, v2];

  % initialize the output cell array.
  X = cell(N, 1);
  for k = 1 : N
    X{k} = zeros(Ka, Kb);
  end

  % loop through the roi boundaries.
  k = 1;
  for idx = 1 : rows(v)
    % get the first dimension roi start and end indices.
    i1 = v(idx, 1);
    i2 = v(idx, 2);

    % get the second dimension roi start and end indices.
    j1 = v(idx, 3);
    j2 = v(idx, 4);

    % get the roi variable count.
    ni = i2 - i1 + 1;
    nj = j2 - j1 + 1;
    nk = ni * nj;

    % extract the segment of interest.
    Xi = Xroi(:, k : k + nk - 1);

    % loop for each element of the cell array.
    for n = 1 : N
      % de-vectorize and insert the current segment.
      Xik = reshape(Xi(n, :), nj, ni);
      X{n}(j1 : j2, i1 : i2) = Xik;
    end

    % shift past the extracted region.
    k += nk;
  end
end

