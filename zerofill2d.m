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
## @anchor{zerofill2d}
## @deftypefn {Function File} {@var{zfid} =} zerofill2d (@var{fid}, @var{k})
## Appends zeros at the end of a two-dimensional free induction decay (FID)
## matrix or cell array by doubling the total length @var{k} times along
## each dimension.
##
## Instead of using this function directly, it is recommended that you use
## @ref{zerofill}.
## @end deftypefn

function zfid = zerofill2d (fid, k)
  % check for proper arguments.
  if (nargin < 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % catch any non-vector type of k that we can accept and transform into the
  % expected data type.
  if (isscalar(k))
    % transform k into a vector.
    k = [k; k];
  elseif (iscell(k) && length(k) == 2)
    % transform k into a vector.
    k = [k{1}; k{2}];
  end

  % make sure k is a two-element vector now.
  if (!isvector(k) || length(k) != 2 || any(k < 0) || any(round(k) != k))
    % invalid type.
    error('zerofill2d: number of zero fills must be non-negative integer(s)');
  end

  % check the type of the input data.
  if (ismatrix(fid))
    % get the old and new column count of the matrix.
    c = columns(fid);
    r = rows(fid);
    nc = c * (2 ^ k(1));
    nr = r * (2 ^ k(2));

    % zero fill the input matrix column-wise.
    sufc = zeros(r, nc - c);
    zfid = [fid, sufc];

    % zero fill the input matrix row-wise.
    sufr = zeros(nr - r, nc);
    zfid = [zfid; sufr];
  elseif (iscell(fid))
    % initialize the output cell array.
    zfid = cell(size(fid));
    
    % apodize each matrix in the cell array.
    for idx = 1 : length(fid)
      % zero fill the matrix.
      zfid{idx} = zerofill2d(fid{idx}, k);
    end
  else
    % invalid data type. throw an exception.
    error('zerofill2d: input data must be a matrix or a cell array');
  end
end

