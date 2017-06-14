## Copyright (C) 2014 University of Nebraska Board of Regents.
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
## @anchor{roi2data}
## @deftypefn {Function File} {@var{Xroi} =} roi2data (@var{X}, @var{ab}, @var{parms}, @var{roi})
## @deftypefnx {Function File} {[@var{Xroi}, @var{abroi}] =} roi2data (@var{X}, @var{ab}, @var{parms}, @var{roi})
## Builds a data matrix from paired one- or two-dimensional spectral data and
## regions of interest by concatenating the spectral data inside each region.
## For one-dimensional data, this is effectively the logical inverse of the
## @ref{rmroi} function. For two-dimensional data, regions are vectorized
## prior to concatenation.
## @end deftypefn

function [Xroi, abroi] = roi2data (X, ab, parms, roi)
  % check the number of input arguments.
  if (nargin != 4 || !any(nargout == [1 : 2]))
    % throw an exception.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are catting one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional cat function.
    [Xroi, abroi] = roi2data1d(X, ab, roi);
  elseif (nd == 2)
    % run a two-dimensional cat function.
    [Xroi, abroi] = roi2data2d(X, ab, roi);
  else
    % throw an exception.
    error('roi2data: input data must be one- or two-dimensional');
  end
end

