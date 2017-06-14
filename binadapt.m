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
## @anchor{binadapt}
## @deftypefn {Function File} {@var{xnew} =} binadapt (@var{X}, @var{ab}, @var{parms}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binadapt (@var{X}, @var{ab}, @var{parms}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binadapt (@var{X}, @var{ab}, @var{parms}, @var{w})
## Adaptively bin a one- or two-dimensional spectrum or spectral dataset in
## @var{X} such that final bins have a width no less than @var{w}. The
## optionally returnable values in @var{abnew} correspond to the new bin
## centers in abscissa units. The final optional return value (@var{widths})
## provides the widths of all the output bins.
##
## In the one-dimensional case, data in @var{X} may either be a column
## vector or a data matrix where each observation is arranged as a row
## in the matrix.
##
## In the two-dimensional case, data in @var{X} may either be a data matrix
## where each direct-dimension slice is along the rows, or a cell array that
## contains multiple matrices, each having direct-dimension slices along
## its rows.
##
## This code is based on the description of `AI-binning' presented in:
##
## @quotation
## De Meyer et. al., `NMR-Based Characterization of Metabolic Alterations in
## Hypertension Using an Adaptive, Intelligent Binning Algorithm',
## Analytical Chemistry, 2008.
## @end quotation
## @end deftypefn

function [xnew, abnew, widths] = binadapt (X, ab, parms, w, R)
  % check for proper arguments.
  if (!any(nargin == [3 : 5]) || !any(nargout == [1 : 3]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % was a width parameter passed?
  if (nargin >= 4 && !isempty(w))
    % check the type.
    if (!isscalar(w) && !isvector(w))
      % invalid type. throw an exception.
      error('binadapt: minimum bin width must be a scalar or a vector');
    end
  else
    % no. use the default.
    w = [];
  end

  % was a resolution parameter passed?
  if (nargin >= 5 && !isempty(R))
    % check the type.
    if (!isscalar(R) || !isreal(R))
      % invalid type. throw an exception.
      error('binadapt: resolution param must be a real scalar');
    end
  else
    % no. use the default.
    R = [];
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are binning one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional binning.
    [xnew, abnew, widths] = binadapt1d(X, ab, w, R);
  elseif (nd == 2)
    % run a two-dimensional binning.
    [xnew, abnew, widths] = binadapt2d(X, ab, w, R);
  else
    % throw an exception.
    error('binadapt: input data must be one- or two-dimensional');
  end
end

