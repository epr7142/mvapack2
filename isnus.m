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
## @anchor{isnus}
## @deftypefn {Function File} {@var{tf} =} isnus (@var{x})
## Returns whether a vector is non-uniformly sampled or not. An abscissa
## vector @var{x} is deemed non-uniformly sampled when any absolute
## difference between consecutive points is greater than or equal to
## twice the minimum difference between all its consecutive points.
## @end deftypefn

function tf = isnus (x)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1 || !isvector(x) || !isreal(x))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % calculate a difference vector to find jumps.
  xdiff = abs([0; x] - [x; 0]);
  xdiff = xdiff(2 : end - 1);
  xjump = sort(find(xdiff >= 2 * min(xdiff)));

  % determine if any jumps were found.
  if (any(xjump))
    % yes! nonuniform.
    tf = true;
  else
    % no. uniform.
    tf = false;
  end
end

