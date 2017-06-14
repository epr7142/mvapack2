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
## @anchor{flatness}
## @deftypefn {Function File} {@var{f} =} flatness (@var{X})
## @deftypefnx {Function File} {@var{f} =} flatness (@var{X}, @var{ab}, @var{roi})
## Calculate the spectral flatness, the ratio of the geometric mean to
## the arithmetic mean of a signal. The optional arguments @var{ab} and
## @var{roi} may be used to specify regions of interest for which to
## calculate flatness, instead of the entire data matrix (default).
##
## If a single vector is provided in @var{X}, one flatness will be returned
## for each ROI. If a data matrix is provided, each observation (row) in the
## matrix will return its own set of flatness values.
## @end deftypefn

function f = flatness (X, ab, roi)
  % check for proper arguments.
  if (!any(nargin == [1 : 3]) || nargout != 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if an ROI matrix was provided.
  if (nargin == 1)
    % no. build a fake abscissa vector and ROI matrix.
    R = [1, max(columns(X), length(X))];
    ab = [R(1) : R(2)];
  elseif (nargin == 3 && isvector(ab) && ismatrix(roi) && columns(roi) == 2)
    % extract the matrix indices from the ROI matrix.
    R = sort(findnearest(ab, roi), 2);
  else
    % invalid argument arrangement. throw an exception.
    print_usage();
  end

  % initialize the output value.
  f = [];

  % determine the type of the input data.
  if (isvector(X))
    % loop through the ROIs.
    for i = 1 : rows(R)
      % extract the subvector of interest.
      Xi = abs(X(R(i, 1) : R(i, 2))) .^ 2;

      % get the current flatness value.
      fi = mean(Xi, 'g') ./ mean(Xi);
      fi = max(min(fi, 1), 0);
      f = [f, fi];
    end
  elseif (ismatrix(X))
    % loop through the ROIs.
    for i = 1 : rows(R)
      % extract the submatrix of interest.
      Xi = abs(X(:, R(i, 1) : R(i, 2))) .^ 2;

      % get the current flatness value.
      fi = mean(Xi, 2, 'g') ./ mean(Xi, 2);
      fi = max(min(fi, 1), 0);
      f = [f, fi];
    end
  else
    % invalid data type. throw an exception.
    error('flatness: data argument must be a vector or matrix');
  end
end

