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
## @anchor{orthspace}
## @deftypefn {Function File} {@var{V} =} orthspace (@var{X}, @var{Y})
## Returns a matrix containing the @var{Y}-orthonormal subspace of
## the data matrix @var{X}. This is used during OPLS modeling, typically.
## @end deftypefn

function V = orthspace (X, Y)
  % check for proper arguments.
  if (nargin != 2 || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % initialize the output matrix.
  V = [];

  % loop through the number of classes (columns of Y).
  for m = 1 : columns(Y)
    % extract the current class matrix column.
    y = Y(:,m);

    % extract the y-orthogonal data from X into v.
    v = (X' * y) ./ (y' * y);
    V = [V, v];
  end
end

