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
## @anchor{realnmr}
## @deftypefn {Function File} {@var{x} =} realnmr (@var{s}, @var{parms})
## Discards imaginaries from all dimensions of an NMR dataset in @var{s} and
## returns only the real spectrum in @var{x}. A parameter structure (or array)
## must be passed as a second argument.
## @end deftypefn

function x = realnmr (s, parms)
  % check if the number of expected arguments was passed.
  if (nargin != 2 || nargout != 1)
    % print the usage statement.
    print_usage();
  end

  % get the dimensions in the dataset.
  nd = nmrdims(s, parms);

  % check whether we are realifying one- or two-dimensional data.
  if (nd == 1)
    % determine the data type of the spectral data.
    if (isvector(s) || ismatrix(s))
      % this is as simple as returning the real component of the data.
      x = real(s);
    else
      % throw an exception.
      error('realnmr: one-dimensional data must be a vector or a matrix');
    end
  elseif (nd == 2)
    % determine the data type of the spectral data.
    if (ismatrix(s))
      % extract the real component of the real states matrix.
      if (iscomplex(s))
        % cut the deck and return the reals.
        x = real(states(s));
      else
        % already real. don't cut the deck again.
        x = s;
      end
    elseif (iscell(s))
      % initialize the output cell array.
      x = cell(length(s), 1);

      % realify each matrix in the cell array.
      for idx = 1 : length(s)
        % realify the matrix.
        if (iscomplex(s{idx}))
          % cut the deck and return the reals.
          x{idx} = real(states(s{idx}));
        else
          % already real. don't cut the deck again.
          x{idx} = s{idx};
        end
      end
    else
      % throw an exception.
      error('realnmr: two-dimensional data must be a matrix or a cell array');
    end
  else
    % throw an exception.
    error('realnmr: input data must be one- or two-dimensional');
  end
end

