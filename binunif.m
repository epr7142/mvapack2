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
## @anchor{binunif}
## @deftypefn {Function File} {@var{xnew} =} binunif (@var{X}, @var{ab}, @var{parms}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binunif (@var{X}, @var{ab}, @var{parms}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binunif (@var{X}, @var{ab}, @var{parms}, @var{w})
## @deftypefnx {Function File} {@var{xnew} =} binunif (@var{X}, @var{ab}, @var{parms}, @var{w})
## Manually bin a one- or two-dimensional spectrum or spectral dataset in
## @var{X} based on uniform bin widths provided in @var{w}. The @var{abnew}
## and @var{widths} values are optionally returnable.
## @end deftypefn

function [xnew, abnew, widths] = binunif (X, ab, parms, w)
  % check for proper arguments.
  if (!any(nargin == [3 : 4]) || !any(nargout == [1 : 3]))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the bin width has been specified.
  if (nargin < 4 || isempty(w))
    % use the default value.
    w = [];
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are binning one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional binning.
    [xnew, abnew, widths] = binunif1d(X, ab, w);
  elseif (nd == 2)
    % run a two-dimensional binning.
    [xnew, abnew, widths] = binunif2d(X, ab, w);
  else
    % throw an exception.
    error('binunif: input data must be one- or two-dimensional');
  end
end

