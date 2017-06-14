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
## @anchor{roiplot}
## @deftypefn {Function File} {} roiplot (@var{X}, @var{ab}, @var{parms}, @var{roi})
## Builds a line plot of one-dimensional spectral data with overlaid regions
## of interest as rectangles, or builds a contour plot of two-dimensional
## spectral data with overlaid regions of interest as rectangles.
## @end deftypefn

function roiplot (X, ab, parms, roi)
  % check the number of input arguments.
  if (nargin != 4)
    % throw an exception.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are plotting one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional plot function.
    roiplot1d(X, ab, roi);
  elseif (nd == 2)
    % run a two-dimensional plot function.
    roiplot2d(X, ab, roi);
  else
    % throw an exception.
    error('roiplot: input data must be one- or two-dimensional');
  end
end

