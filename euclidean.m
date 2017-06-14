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
## @anchor{euclidean}
## @deftypefn {Function File} {@var{retval} =} euclidean (@var{X}, @var{y})
## Return the squared Euclidean distance between the means of the multivariate
## samples @var{x} and @var{y}, which must have the same number of
## components (columns), but may have a different number of observations
## (rows).
## @end deftypefn

function retval = euclidean (x, y)
  % check if the number of arguments is correct.
  if (nargin != 2 || nargout != 1)
    % nope. throw an exception.
    error('euclidean: invalid argument count');
  end

  % check that the input matrices contain the proper data type.
  if (!(isnumeric(x) || islogical(x)) || !(isnumeric(y) || islogical(y)))
    % nope. throw an exception.
    error('euclidean: X and Y must be numeric matrices or vectors');
  end

  % check if the input matrices are two-dimensional.
  if (ndims(x) != 2 || ndims(y) != 2)
    % nope. throw an exception.
    error('euclidean: X and Y must be 2-D matrices or vectors');
  end

  % get the dimensions of the input matrices.
  [xr, xc] = size(x);
  [yr, yc] = size(y);

  % check if the dimensionality of the data matches.
  if (xc != yc)
    % mismatch. throw an exception.
    error ('euclidean: X and Y must have the same number of columns');
  end

  % transform the first input matrix to double-precision, if necessary.
  if (isinteger(x))
    x = double(x);
  end

  % calculate the means of the data matrices.
  xm = mean(x);
  ym = mean(y);

  % subtract the means from the matrices.
  x = bsxfun (@minus, x, xm);
  y = bsxfun (@minus, y, ym);

  % return the result.
  retval = (xm - ym) * (xm - ym)';
end

