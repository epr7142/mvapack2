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
## @anchor{roifun}
## @deftypefn {Function File} {@var{Y} =} roifun (@var{X}, @var{ab}, @var{parms}, @var{roi}, @var{func})
## Executes a function handle @var{func} within each region of interest in a
## one- or two-dimensional spectral dataset and returns the resulting values
## in @var{Y}. The value produced by @var{func} must be scalar, or this
## function will not execute it.
## @end deftypefn

function Y = roifun (X, ab, parms, roi, func)
  % check the number of input arguments.
  if (nargin != 5 || nargout != 1)
    % throw an exception.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(X, parms);

  % check whether we are mapping one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional map function.
    Y = roifun1d(X, ab, roi, func);
  elseif (nd == 2)
    % run a two-dimensional map function.
    Y = roifun2d(X, ab, roi, func);
  else
    % throw an exception.
    error('roifun: input data must be one- or two-dimensional');
  end
end

