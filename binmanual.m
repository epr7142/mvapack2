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
## @anchor{binmanual}
## @deftypefn {Function File} {@var{xnew} =} binmanual (@var{X}, @var{ab}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binmanual (@var{X}, @var{ab}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binmanual (@var{X}, @var{ab}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {@var{xnew} =} binmanual (@var{X}, @var{ab}, @var{parms}, @var{centers}, @var{widths})
## Manually bin a one-dimensional spectrum or spectral dataset in @var{X} based
## on regions of interest provided in @var{roi}, or centers and widths provided
## in @var{centers} and @var{widths}. If regions of interest are used to bin,
## the @var{abnew} and @var{widths} values are optionally returnable.
## @end deftypefn

function [xnew, abnew, widths] = binmanual (X, ab, parms, a, b)
  % check for proper arguments.
  if (!any(nargin == [4 : 5]) || !any(nargout == [1 : 3]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check which argument structure was provided.
  if (nargin == 4)
    % get the roi matrix.
    roi = a;
  elseif (nargin == 5)
    % get the roi matrix.
    roi = bin2roi(a, b);
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are binning one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional binning.
    [xnew, abnew, widths] = binmanual1d(X, ab, roi);
  elseif (nd == 2)
    % run a two-dimensional binning.
    [xnew, abnew, widths] = binmanual2d(X, ab, roi);
  else
    % throw an exception.
    error('binmanual: input data must be one- or two-dimensional');
  end
end

