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
## @anchor{decompose}
## @deftypefn {Function File} {@var{T} =} decompose (@var{f}, @var{t}, @var{parms})
## @deftypefnx {Function File} {@var{T} =} decompose (@var{f}, @var{t}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {@var{T} =} decompose (@var{f}, @var{t}, @var{parms}, @var{roi}, @var{minbw})
## @deftypefnx {Function File} {@var{T} =} decompose (@var{f}, @var{t}, @var{parms}, @var{roi}, @var{minbw}, @var{fitopts})
## Performs Complete Reduction to Amplitude and Frequency Table (CRAFT)
## analysis of a time-domain NMR free induction decay (FID) in @var{f},
## with a time abscissa in @var{t} and parameters in @var{parms}. If
## multiple decays are provided in @var{f}, then a joined table
## reflecting data from all decays will be returned.
##
## In the one-dimensional case, data in @var{f} may either be a column
## vector or a data matrix where each observation is arranged as a row
## in the matrix.
##
## In the two-dimensional case, data in @var{f} may either be a data matrix
## where each direct-dimension slice is along the rows, or a cell array that
## contains multiple matrices, each having direct-dimension slices along
## its rows.
##
## The optional argument @var{roi} may be passed to specify either a function
## handle for automated region of interest (ROI) selection, or a matrix of
## manually defined regions of interest. Each manually defined region should be
## a two-element row in @var{roi} containing the lower and upper frequency
## values (in hertz, @xref{nmrft}). Two-dimensional data should have ROI rows
## with four numbers (max and min for each dimension).
## In the case of a function handle, @var{roi} must be a function
## defined as follows:
##
## @code{function roi = roi_function (s, ab, parms, wmin)
##   ...
## end}
##
## The optional argument @var{minbw} may be passed to specify a minimum
## bandwidth (in Hertz) of the automatically selected ROIs. The default value
## of @var{minbw} is one one-hundredth of the total spectral width.
##
## The CRAFT algorithm is implemented according to information reported in:
##
## @quotation
## Krishnamurthy K., `CRAFT (complete reduction to amplitude frequency table) -
## robust and time-efficient Bayesian approach for quantitative mixture
## analysis by NMR', Magnetic Resonance in Chemistry, 2013.
## @end quotation
## @end deftypefn

function T = decompose (f, t, parms, roi, minbw, fitopts)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [3 : 6]) || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % check if an ROI argument was provided.
  if (nargin < 4 || isempty(roi))
    % use the default automatic ROI generation routine.
    roi = [];
  end

  % check if a minimum bandwidth was provided.
  if (nargin < 5 || isempty(minbw))
    % use the default minimum bandwidth.
    minbw = [];
  end

  % check if an options structure was passed.
  if (nargin < 6 || isempty(fitopts))
    % pass an empty options structure to the fitting routine.
    fitopts = [];
  end

  % get the dimensions in the dataset.
  nd = nmrdims(f, parms);

  % check whether we are decomposing one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional decomposition.
    T = decompose1d(f, t, parms, roi, minbw, fitopts);
  elseif (nd == 2)
    % run a two-dimensional decomposition.
    T = decompose2d(f, t, parms, roi, minbw, fitopts);
  else
    % throw an exception.
    error('decompose: input data must be one- or two-dimensional');
  end
end

