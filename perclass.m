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
## @anchor{perclass}
## @deftypefn {Function File} {@var{Xnew} =} perclass (@var{fn}, @var{X}, @var{Y})
## Performs the operation defined in the function handle @var{fn} on @var{X}
## in a class-dependent manner. In other words, instead of treating the entire
## dataset as one as would occur when applying @var{fn} normally, @var{fn} is
## applied to each class @emph{individually} and the results are reassembled
## into the matrix @var{Xnew}.
##
## It is important to note that the function @var{fn} must not modify either
## the number of rows or columns of the data matrix @var{X}. Finally, the
## function handle @var{fn} must be of the following form:
## @example
## function Xnew = fn (X), ... end
## @end example
## @end deftypefn

function Xnew = perclass (fn, X, Y)
  % check for proper arguments.
  if (nargin != 3 || nargout != 1 || !is_function_handle(fn) || ...
      !ismatrix(X) || !ismatrix(Y))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check that the observation counts match.
  if (rows(X) != rows(Y))
    % throw an exception.
    error('perclass: observation counts of X and Y do not match');
  end

  % initialize the output matrix.
  Xnew = zeros(size(X));

  % loop through the classes in the class matrix.
  for m = 1 : columns(Y)
    % get the indices of the current class.
    idx = classidx(Y, m);

    % apply the function to the observations of the current class.
    Xnew(idx,:) = fn(X(idx,:));
  end
end

