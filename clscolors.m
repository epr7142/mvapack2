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
## @anchor{clscolors}
## @deftypefn {Function File} {@var{colors} =} clscolors (@var{Y})
## Uses the @var{Y} matrix created by @ref{classes} to build different
## colors for each class.
## @end deftypefn

function [colors, clrmap] = clscolors (Y)
  % check the type of arguments.
  if (nargin != 1 || !any(nargout == [1 : 2]) || !ismatrix(Y))
    % invalid arguments. throw an exception.
    print_usage();
  end

  % extract the number of classes from Y.
  M = columns(Y);

  % build a colormap from an HSV colorwheel.
  cmap = hsv(M + 1);
  cmap = cmap(1 : M, :);

  % multiply the class matrix by the colormap to produce colors for all
  % observations. this will also work with continuous class membership,
  % i.e. PLS regression instead of discriminant analysis.
  colors = Y * cmap;

  % see if the user requested the original map too.
  if (nargout >= 2)
    clrmap = cmap;
  end
end

