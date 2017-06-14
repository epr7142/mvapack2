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
## @anchor{roibin}
## @deftypefn {Function File} {@var{roi} =} roibin (@var{s}, @var{ab}, @var{parms}, @var{wmin})
## Generate regions of interest (ROIs) from a spectrum or set of spectra in
## @var{s} with corresponding abscissa in @var{ab}. The minimum viable ROI
## width is specified in @var{wmin}.
##
## This function generates ROIs by first using @ref{binadapt} to bin the
## variables of a spectrum or a set of spectra, and converts the bins into
## ROIs.
## @end deftypefn

function roi = roibin (s, ab, parms, wmin)
  % check for proper arguments.
  if (nargin != 4 || nargout != 1 || ...
      !isvector(ab) || !(isscalar(wmin) || isvector(wmin)) || ...
      (!ismatrix(s) && !isvector(s) && !iscell(s)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % grab the real portion of the spectrum.
  if (iscomplex(s))
    % only if complex values exist in the input.
    x = realnmr(s, parms);
  else
    % if the input data is real, assume we've removed reals already.
    x = s;
  end

  % call the adaptive binning routine.
  [snew, abnew, widths] = binadapt(x, ab, parms, wmin);

  % translate the bins into regions of interest.
  roi = bin2roi(abnew, widths);

  % get the number of dimensions in the data.
  nd = nmrdims(s, parms);

  % estimate the baseline noise level of the data.
  [mu, sigma] = estnoise(x, parms);
  nf = mu + 6 .* sigma;

  % initialize a list of 'keepers'. these regions will be deemed acceptable
  % based on signal to noise criteria.
  keepers = zeros(rows(roi), 1);

  % loop through the regions of interest.
  for idx = 1 : rows(roi)
    % get the min and max values of the ROI for each dimension.
    Fmin = [];
    Fmax = [];
    for d = 0 : nd - 1
      Fmin = [Fmin; roi(idx, 1 + 2 * d)];
      Fmax = [Fmax; roi(idx, 2 + 2 * d)];
    end

    % extract the subspectrum corresponding to the region of interest.
    [xsub, absub] = subspect(x, ab, parms, Fmin, Fmax);

    % grab the spectral maximum based on the dimensionality.
    if (iscell(xsub))
      % grab the maximum from all the spectra.
      xmax = max(cellfun(@(M) max(vec(M)), xsub));
    elseif (isvector(xsub) || ismatrix(xsub))
      % just grab the maximum value.
      xmax = max(vec(xsub));
    end

    % see if the maximum value in the bin exceeds the noise floor.
    if (xmax > median(nf))
      % label this region as a keeper.
      keepers(idx) = 1;
    end
  end

  % just keep the acceptable regions of interest.
  roi = roi(find(keepers),:);
end

