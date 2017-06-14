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
## @anchor{bin2roi}
## @deftypefn {Function File} {@var{roi} =} bin2roi (@var{abnew}, @var{widths})
## Translate bin centers and widths in @var{abnew} and @var{widths},
## respectively, into regions of interest in @var{roi}. The resulting
## regions may then be plotted using @ref{roiplot}.
## @end deftypefn

function roi = bin2roi (abnew, widths)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check the type of the bin widths argument.
  if (isscalar(abnew) && isscalar(widths))
    % calculate the single roi.
    roi = [abnew, abnew] + 0.5 .* widths .* [-1, 1];
  elseif (isvector(abnew) && isvector(widths))
    % ensure the bin vectors have equal length.
    if (length(abnew) != length(widths))
      % length mismatch. throw an exception.
      error('bin2roi: bin vector lengths must match');
    end

    % get the vector lengths.
    n = length(abnew);

    % ensure the bins are column vectors.
    abnew = reshape(abnew, n, 1);
    widths = reshape(widths, n, 1);

    % calculate the roi matrix.
    roi = [abnew, abnew] + 0.5 .* (widths * [-1, 1]);
  elseif (ismatrix(abnew) && ismatrix(widths))
    % ensure the bin matrices have equal row counts.
    if (rows(abnew) != rows(widths))
      % row count mismatch. throw an exception.
      error('bin2roi: bin matrix row counts must match');
    end

    % get the matrix row counts.
    n = rows(abnew);

    % initialize the roi matrix.
    roi = [];

    % loop for every dimension in the bin matrices.
    for idx = 1 : columns(abnew)
      % calculate the current dimension's roi matrix.
      roi = [roi, abnew(:,idx) * [1, 1] + 0.5 .* (widths(:,idx) * [-1, 1])];
    end
  else
    % invalid type. throw an exception.
    error('bin2roi: bins must be a scalar, vector or matrix arguments');
  end
end

