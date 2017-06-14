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
## @anchor{zerofill1d}
## @deftypefn {Function File} {@var{zfid} =} zerofill1d (@var{fid}, @var{k})
## Appends zeros at the end of a one-dimensional free induction decay (FID)
## vector or matrix by doubling the total length @var{k} times.
##
## Instead of using this function directly, it is recommended that you use
## @ref{zerofill}.
## @end deftypefn

function zfid = zerofill1d (fid, k)
  % check for proper arguments.
  if (nargin < 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the type of k.
  if (!isscalar(k) || k < 0 || round(k) != k)
    % invalid type. throw an exception.
    error('zerofill1d: number of zero fills must be a non-negative integer');
  end

  % check the type of the input data.
  if (isvector(fid))
    % ensure the vector is a column vector.
    fid = reshape(fid, length(fid), 1);

    % get the old and new length of the vector.
    r = length(fid);
    n = r * (2 ^ k);

    % zero fill the input vector.
    suffix = zeros(n - r, 1);
    zfid = [fid; suffix];
  elseif (ismatrix(fid))
    % get the old and new column count of the matrix.
    c = columns(fid);
    r = rows(fid);
    n = c * (2 ^ k);

    % zero fill the input matrix.
    suffix = zeros(r, n - c);
    zfid = [fid, suffix];
  else
    % invalid data type. throw an exception.
    error('zerofill1d: input data must be a vector or a matrix');
  end
end

