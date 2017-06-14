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
## @anchor{roipeak}
## @deftypefn {Function File} {@var{roi} =} roipeak (@var{s}, @var{ab}, @var{parms}, @var{wmin})
## Generate regions of interest (ROIs) from a spectrum or set of spectra in
## @var{s} with corresponding abscissa in @var{ab}. The minimum viable ROI
## width is specified in @var{wmin}.
##
## This function generates ROIs by first using @ref{peakpick} to pick the
## peaks of a spectrum or a set of spectra, creates ROIs of width @var{wmin}
## for each peak, and then joins all overlapped ROIs together until no two
## ROIs overlap.
## @end deftypefn

function roi = roipeak (s, ab, parms, wmin)
  % check for proper arguments.
  if (nargin != 4 || nargout != 1 || ...
      !isvector(ab) || !isscalar(wmin) || ...
      (!ismatrix(s) && !isvector(s)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % define a matrix to contain all picked peaks.
  P = [];

  % determine what type of spectral data was passed.
  if (isvector(s))
    % run a single peak pick invocation.
    P = peakpick(s, ab);
  elseif (ismatrix(s))
    % loop through the observations.
    for i = 1 : rows(s)
      % run a peak pick invocation for the current observation.
      Pi = peakpick(s(i,:)', ab);

      % append the picked peaks to the current list.
      P = [P; Pi];
    end
  end

  % generate the initial ROI list.
  roi = P(:,1) * [1, 1] + ones(rows(P), 1) * (0.5 .* wmin .* [-1, 1]);

  % loop through the ROIs.
  idx = 1;
  while (idx < rows(roi))
    % check if the current ROI overlaps the next.
    if (roi(idx, 2) >= roi(idx + 1, 1))
      % overlap. merge and delete the next ROI.
      roi(idx, 2) = roi(idx + 1, 2);
      roi(idx + 1, :) = [];

      % re-sort the ROI list.
      [jnk, newidx] = sort(roi(:,1));
      roi = roi(newidx, :);
    else
      % no overlap. move to the next ROI.
      idx++;
    end
  end
end

