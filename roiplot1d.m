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
## @anchor{roiplot1d}
## @deftypefn {Function File} {} roiplot1d (@var{roi})
## @deftypefnx {Function File} {} roiplot1d (@var{X}, @var{ab}, @var{roi})
## Overlays regions of interest (as lines) on an existing line plot of
## frequency-domain spectral data, or builds a line plot of data with
## overlaid regions of interest as rectangles.
## @end deftypefn

function roiplot1d (a, b, c)
  % check the number of input arguments.
  if (nargin == 1 && ismatrix(a) && columns(a) == 2)
    % we are overlaying ROI lines.
    havex = false;
    roi = a;
  elseif (nargin == 3 && ismatrix(a) && isvector(b) && ...
          ismatrix(c) && columns(c) == 2)
    % we are overlaying ROI rectangles.
    havex = true;
    X = real(a);
    ab = b;
    roi = c;
  else
    % unknown. throw an exception.
    print_usage();
  end

  % do we have data to plot first?
  if (havex == true)
    % plot the rows of X.
    figure();
    hold on;
    plot(ab, X);
    hold off;
  end

  % get the current axis handle.
  hax = gca();
  sax = get(hax);

  % determine how to grab y-extents for the rectangles.
  if (havex == true)
    % vector or matrix data?
    if (isvector(X))
      % use the extents of the vector, plus some slack.
      yext = [X, X] + ones(size(X)) * (0.01 .* range(X) .* [-1, 1]);
    elseif (ismatrix(X))
      % use the extents of the plotted data matrix.
      yext = [min(X); max(X)]';
      yextrange = max(max(X)) - min(min(X));
      yext += ones(rows(yext), 1) * (0.01 .* yextrange .* [-1, 1]);
    end
  else
    % use the extents of the plot.
    yext = sax.ylim + 0.01 .* [range(sax.ylim), range(sax.ylim)];
  end

  % loop through the regions of interest.
  for i = 1 : rows(roi)
    % get the x-positions and width of the region.
    L = min(roi(i,:));
    R = max(roi(i,:));
    W = range(roi(i,:));

    % determine how to get the y-positions and height.
    if (havex == true)
      % get the abscissa indices of the x-positions.
      abL = min(findnearest(ab, [L, R]));
      abR = max(findnearest(ab, [L, R]));

      % use the (frequency-dependent) data extents.
      B = min(yext([abL : abR], 1));
      T = max(yext([abL : abR], 2));
      H = T - B;
    else
      % use the plot extents.
      B = yext(1);
      T = yext(2);
      H = T - B;
    end

    % plot the rectangle.
    rectangle(hax, 'Position', [L, B, W, H]);
  end
end

