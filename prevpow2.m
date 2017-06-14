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
## @anchor{prevpow2}
## @deftypefn {Function File} {} prevpow2 (@var{x})
## If @var{x} is a scalar, return the last integer N such that
## @code{2^n <= abs(x)}.
##
## If @var{x} is a vector, return @code{prevpow2(length(x))}.
## @end deftypefn

function y = prevpow2 (x)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the input type.
  if (isscalar(x))
    % return the value.
    y = floor(log2(abs(x)));
  elseif (isvector(x))
    % return based on the length of the vector.
    y = prevpow2(length(x));
  else
    % invalid type. output an error.
    error('prevpow2: invalid input data type');
  end
end

