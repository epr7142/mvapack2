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
## @anchor{loadascii}
## @deftypefn {Function File} {@var{x} =} loadascii (@var{filename})
## @deftypefnx {Function File} {[@var{x}, @var{ab}] =} loadascii (@var{filename})
## @deftypefnx {Function File} {@var{X} =} loadascii (@var{filenames})
## @deftypefnx {Function File} {[@var{X}, @var{ab}] =} loadascii (@var{filenames})
## Loads one or more ASCII files, each which contains a two-column,
## space-delimited format. The first column is expected to be an abscissa
## and the second column is expected to be the data.
## @end deftypefn

function [X, ab] = loadascii (filenames)
  % check for the minimum argument count.
  if (nargin != 1 || !any(nargout == [1 : 2]))
    % this won't work. throw an error.
    print_usage();
  end

  % check the type of the file input argument.
  if (ischar(filenames))
    % turn the single string into a cell array...
    filenames = {filenames};
  elseif (iscellstr(filenames))
    % make sure at least one file was passed.
    if (length(filenames) == 0)
      % no. throw an exception.
      error('loadascii: zero input files supplied');
    end
  else
    % invalid input type. throw an exception.
    error('loadascii: filenames must be a string or a string cell array');
  end

  % initialize the output data matrix and observation count.
  n = length(filenames);
  ax = [];
  X = [];

  % allocate a temporary cell array.
  C = cell(n, 1);

  % loop through the filenames.
  for i = 1 : n
    % try to load in the text file.
    try
      C{i} = load(filenames{i});
    catch
      error('loadascii: failed to load file "%s"', filenames{i});
    end

    % ensure the number of columns is correct.
    if (columns(C{i}) != 2)
      % throw an exception.
      error('loadascii: expected two-column data in "%s"', filenames{i});
    end

    % extract columns from the current matrix.
    axi = C{i}(:,1);
    xi = C{i}(:,2);

    % determine which observation this is.
    if (i == 1)
      % the first. store the axis boundaries.
      abmin = min(axi);
      abmax = max(axi);

      % construct a uniform axis from the boundaries.
      ax = [abmin : min(diff(axi)) : abmax];
    else
      % the second or later. check the axis.
      if (min(axi) != abmin || max(axi) != abmax)
        % throw an exception.
        error('loadascii: abscissa of "%s" does not match', filenames{i});
      end
    end

    % interpolate the values of the observation and store it.
    xi = interp1(axi, xi, ax, 'linear');
    X = [X; reshape(xi, 1, length(xi))];
  end

  % see if the output data should be reshaped.
  if (n == 1)
    % yes. reshape the data into a column vector.
    X = reshape(X, length(xi), 1);
  end

  % see if the user requested the axis to be returned.
  if (nargout >= 2)
    % return the axis.
    ab = reshape(ax, length(ax), 1);
  end
end

