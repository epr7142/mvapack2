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
## @anchor{findnearest}
## @deftypefn {Function File} {@var{idx} =} findnearest (@var{x}, @var{a})
## Finds the nearest index in a vector @var{x} to a value @var{a}. This is
## similar to the @code{find} command, but returns the nearest match if no
## entries are exact.
## @end deftypefn

function idx = findnearest (x, a)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1 || !isvector(x))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the type of the search value.
  if (isscalar(a))
    % return the nearest index.
    idx = find(min(abs(x - a)) == abs(x - a));
  elseif (ismatrix(a))
    % return a matrix of all the found indices.
    idx = arrayfun(@(v) findnearest(x, v), a);
  else
    % invalid type.
    error('findnearest: search value "a" must be a scalar or matrix');
  end
end

