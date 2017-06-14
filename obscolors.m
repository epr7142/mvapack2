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
## @anchor{obscolors}
## @deftypefn {Function File} {@var{colors} =} obscolors (@var{XY})
## Uses the @var{XY} matrix (either @var{X} or @var{Y}) to build
## different colors for each observation, rainbow-style.
## @end deftypefn

function colors = obscolors (XY)
  % check the type of arguments.
  if (nargin != 1 || nargout != 1 || !ismatrix(XY))
    % invalid arguments. throw an exception.
    print_usage();
  end

  % get the number of observations.
  N = rows(XY);

  % calculate the color map from the HSV colorwheel.
  colors = hsv(N + 1);
  colors = colors(1 : N, :);
end

