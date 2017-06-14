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
## @anchor{roiplot2d}
## @deftypefn {Function File} {} roiplot2d (@var{roi})
## @deftypefnx {Function File} {} roiplot2d (@var{X}, @var{ab}, @var{roi})
## Overlays regions of interest (as rectangles) on an existing contour plot of
## frequency-domain spectral data, or builds a contour plot of data with
## overlaid regions of interest as rectangles.
## @end deftypefn

function roiplot2d (a, b, c)
  % check the number of input arguments.
  if (nargin == 1 && ismatrix(a) && columns(a) == 4)
    % we are overlaying ROI lines.
    havex = false;
    roi = a;
  elseif (nargin == 3 && (ismatrix(a) || iscell(a)) && iscell(b) && ...
          ismatrix(c) && columns(c) == 4)
    % we are overlaying ROI rectangles.
    havex = true;
    X = a;
    ab = b;
    roi = c;
  else
    % unknown. throw an exception.
    print_usage();
  end

  % do we have data to plot first?
  if (havex == true)
    % check the data matrix.
    if (iscomplex(X))
      % throw an exception.
      error('roiplot2d: input data matrix must not be complex');
    end

    % build a mesh grid of the abscissa values.
    [gx, gy] = meshgrid(ab{1}, ab{2});

    % start plotting.
    figure();
    hold on;

    % determine how to plot X.
    if (iscell(X))
      % build a matrix of the average dataset intensity.
      Z = zeros(size(X{1}));
      for n = 1 : length(X)
        Z += X{1} ./ n;
      end

      % plot the contours of the average matrix.
      contour(gx, gy, Z);
    else
      % plot the contours of the matrix.
      contour(gx, gy, X);
    end

    % finish plotting.
    hold off;
  end

  % get the current axis handle.
  hax = gca();
  sax = get(hax);

  % loop through the regions of interest.
  for i = 1 : rows(roi)
    % get the x-positions and width of the region.
    L = min(roi(i, [1, 2]));
    B = min(roi(i, [3, 4]));
    W = range(roi(i, [1, 2]));
    H = range(roi(i, [3, 4]));

    % plot the rectangle.
    rectangle(hax, 'Position', [L, B, W, H]);
  end
end

