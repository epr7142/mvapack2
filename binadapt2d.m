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
## @anchor{binadapt2d}
## @deftypefn {Function File} {@var{xnew} =} binadapt2d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binadapt2d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binadapt2d (@var{X}, @var{ab}, @var{w})
## Adaptively bin a two-dimensional spectrum or spectral dataset.
##
## It is highly recommended that you use @ref{binadapt} instead of this
## function directly.
## @end deftypefn

function [xnew, abnew, widths] = binadapt2d (X, ab, w, R)
  % check for proper arguments.
  if (!any(nargin == [2 : 4]) || nargout != 3 || ...
      !(ismatrix(X) || iscell(X)) || !iscell(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a third argument was supplied.
  if (nargin < 3 || isempty(w))
    % no. use the default initial bin width of 0.025 abscissa units.
    w = [0.025; 2.0];
  else
    % check the type.
    if (!isvector(w) || rows(w) != 2)
      % invalid type.
      error('binadapt2d: width argument must be a two-element column vector');
    end
  end

  % check if a fourth argument was supplied.
  if (nargin < 4 || isempty(R))
    % no. use the default initial resolution parameter.
    R = 1.0;
  else
    % check the type.
    if (!isscalar(R) || !isreal(R) || R < 0.0 || R > 2.0)
      % invalid type.
      error('binadapt2d: resolution param must be a real scalar in [0,2]');
    end
  end

  % check if the data is a single spectrum.
  if (ismatrix(X))
    % initialize a single-element cell array.
    C = cell(1, 1);
    C{1} = X ./ max(vec(X));
  elseif (iscell(X))
    % use the cell array directly.
    C = cellfun(@(x) x ./ max(vec(x)), X, 'UniformOutput', false);
  else
    % invalid type.
    error('binadapt2d: input data must be a matrix or a cell array');
  end

  % initialize output variables.
  v = [];
  abnew = [];
  widths = [];

  % calculate the minimum bin size in points.
  abi = cellfun(@(x) length(x) / range(x), ab);
  kmin = max(floor(w .* abi), [8; 4]);

  % run the adaptive intelligent binning subroutine.
  v = __binadapt2d(C, round(kmin(2)), round(kmin(1)), R) + 1;

  % initialize the output bin matrix.
  xnew = zeros(length(C), rows(v));

  % loop through the bin boundaries.
  for idx = 1 : rows(v)
    % get the first-dimension bin start and end indices.
    j1 = v(idx, 3);
    j2 = v(idx, 4);

    % get the second-dimension bin start and end indices.
    i1 = v(idx, 1);
    i2 = v(idx, 2);

    % build the abscissas of the current bin.
    ab1 = ab{1}(j1 : j2);
    ab2 = ab{2}(i1 : i2);

    % build the next abscissa center, width and index.
    abnew = [abnew; median(ab1), median(ab2)];
    widths = [widths; abs(range(ab1)), abs(range(ab2))];

    % loop for each element of the cell array.
    for k = 1 : length(C)
      % get the bin segment from the current matrix.
      Ck = C{k};
      Xi = Ck(i1 : i2, j1 : j2);

      % bin the segment into its row.
      xnew(k,idx) = abs(trapz(ab2, trapz(ab1, Xi, 2)));
    end
  end

  % see if the output bin matrix is really a vector.
  if (length(C) == 1)
    % yes. reshape it into column form.
    xnew = reshape(xnew, length(xnew), 1);
  end
end

