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
## @anchor{binunif2d}
## @deftypefn {Function File} {@var{xnew} =} binunif2d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binunif2d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binunif2d (@var{X}, @var{ab}, @var{w})
## Uniformly bin a two-dimensional spectrum or spectral dataset.
##
## It is highly recommended that you use @ref{binunif} instead of this
## function directly. Better yet, stop uniformly binning your datasets
## and use @ref{binadapt}.
## @end deftypefn

function [xnew, abnew, widths] = binunif2d (X, ab, w)
  % check for proper arguments.
  if (!any(nargin == [2 : 3]) || nargout != 3 || ...
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
      error('binunif2d: width argument must be a two-element column vector');
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
    error('binunif2d: input data must be a matrix or a cell array');
  end

  % get the size of the data.
  [Ka, Kb] = size(C{1});

  % initialize output variables.
  xnew = [];
  abnew = [];
  widths = [];

  % calculate the minimum bin size in points.
  abi = cellfun(@(x) length(x) / range(x), ab);
  kbin = max(floor(w .* abi), [8; 4]);

  % loop through all the bins to be made.
  ka = 1;
  while (ka < Ka)
    % move to the next first-dimension bin index.
    ka1 = ka;
    ka2 = min(Ka, ka + kbin(2));

    % skip zero-width bins.
    if (ka1 == ka2)
      ka = ka2;
      continue
    end

    % loop through the second dimension of bins.
    kb = 1;
    while (kb < Kb)
      % move to the next second-dimension bin index.
      kb1 = kb;
      kb2 = min(Kb, kb + kbin(1));

      % skip zero-width bins.
      if (kb1 == kb2)
        kb = kb2;
        continue
      end

      % build the abscissas of the current bin.
      ab1 = ab{1}(kb1 : kb2);
      ab2 = ab{2}(ka1 : ka2);

      % compute the bin center and width.
      abnew = [abnew; median(ab1), median(ab2)];
      widths = [widths; abs(range(ab1)), abs(range(ab2))];

      % loop over the elements of the cell array.
      for n = 1 : length(C)
        % get the bin segment from the current matrix.
        Cn = C{n};
        Xi = Cn(ka1 : ka2, kb1 : kb2);

        % bin the segment into its row.
        xnew = [xnew, abs(trapz(ab2, trapz(ab1, Xi, 2)))];
      end

      % increment the inner loop counter.
      kb = kb2;
    end

    % increment the outer loop counter.
    ka = ka2;
  end

  % see if the output bin matrix is really a vector.
  if (length(C) == 1)
    % yes. reshape it into column form.
    xnew = reshape(xnew, length(xnew), 1);
  end
end

