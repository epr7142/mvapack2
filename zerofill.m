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
## @anchor{zerofill}
## @deftypefn {Function File} {@var{zfid} =} zerofill (@var{fid}, @var{parms})
## @deftypefnx {Function File} {@var{zfid} =} zerofill (@var{fid}, @var{parms}, @var{k})
## Appends zeros at the end of a free induction decay (FID) vector, matrix or
## cell array by doubling the total length @var{k} times. If @var{k} is not
## supplied, a default number of one zero fill is performed.
##
## For one-dimensional data, @var{k} must be a scalar value. In the case of
## two-dimensional data, @var{k} may be a scalar or a vector.
## @end deftypefn

function zfid = zerofill (fid, parms, k)
  % check for proper arguments.
  if (nargin < 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % see if the number of zero fills was passed.
  if (nargin < 3 || isempty(k))
    % no. default to one zero fill.
    k = 1;
  end

  % get the dimensions in the dataset.
  nd = nmrdims(fid, parms);

  % check whether we are zero-filling one- or two-dimensional data.
  if (nd == 1)
    % run a one-dimensional zero-fill.
    zfid = zerofill1d(fid, k);
  elseif (nd == 2)
    % run a two-dimensional zero-fill.
    zfid = zerofill2d(fid, k);
  else
    % throw an exception.
    error('zerofill: input data must be one- or two-dimensional');
  end
end

