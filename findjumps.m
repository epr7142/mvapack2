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
## @anchor{findjumps}
## @deftypefn {Function File} {@var{k} =} findjumps (@var{x})
## In a vector @var{x} that is expected to contain mostly uniformly spaced
## data points, find the indices where the values jump across larger than
## expected regions.
##
## This function is highly accepting of sampling jitter, as it only registers
## a jump when the difference between two consecutive points exceeds twice
## the standard uniform spacing.
## @end deftypefn

function k = findjumps (x)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1 || !isvector(x))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % calculate a difference vector to find jumps.
  xdiff = abs([0; x] - [x; 0]);
  xdiff = xdiff(2 : end - 1);
  xjump = sort(find(xdiff > 2 * xdiff(1)));
  k = xjump + 1;
end

