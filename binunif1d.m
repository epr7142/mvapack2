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
## @anchor{binunif1d}
## @deftypefn {Function File} {@var{xnew} =} binunif1d (@var{X}, @var{ab})
## @deftypefnx {Function File} {@var{xnew} =} binunif1d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binunif1d (@var{X}, @var{ab})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binunif1d (@var{X}, @var{ab}, @var{w})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binunif1d (@var{X}, @var{ab})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binunif1d (@var{X}, @var{ab}, @var{w})
## Uniformly bin a one-dimensional spectrum or spectral dataset in @var{X}
## such that final bins have a width no greater than @var{w}. The optionally
## returnable values in @var{abnew} correspond to the new bin centers in
## abscissa units.
##
## The @var{w} value is optional and has a default value of @var{0.025}.
##
## This function is known to create bin boundaries that result in correlated
## output variables. Unless uniform bins are a requirement of the task at
## hand, @ref{binoptim} is the recommended binning method.
## @end deftypefn

function [xnew, abnew, widths] = binunif1d (X, ab, w)
  % check for proper arguments.
  if (!any(nargin == [2 : 3]) || !any(nargout == [1 : 3]) || ...
      !ismatrix(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a third argument was supplied.
  if (nargin < 3 || isempty(w))
    % no. use the default initial bin width of 0.025 abscissa units.
    w = 0.025;
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('binunif1d: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % initialize variables used during binning.
  kbin = floor(K / floor(range(ab) / w));

  % initialize output variables.
  xnew = [];
  abnew = [];

  % find the jump locations in the abscissa.
  abjump = findjumps(ab) - 1;

  % loop through all the bins to be made.
  k = 1;
  while (k < K)
    % move to the next bin index.
    k1 = k;

    % check if any jumps exist.
    if (length(abjump))
      % increment to avoid the jumps if they are in this bin.
      k2 = min([K; k + kbin; abjump(1)]);

      % remove the avoided jump.
      if (k2 == abjump(1))
        abjump = abjump(2 : end);
      end
    else
      % no jumps remain. increment as usual.
      k2 = min([K; k + kbin]);
    end

    % skip zero-width bins.
    if (k1 == k2)
      k = k2 + 1;
      continue
    end

    % build the segments over which to calculate the next bin.
    Xcut = X(:,k1 : k2);
    abcut = ab(k1 : k2);

    % build the next output bin and abscissa center.
    xnew = [xnew, abs(trapz(abcut, Xcut, 2))];
    abnew = [abnew; median(abcut)];

    % increment the loop counter.
    k = k2 + 1;
  end

  % see if a third output was requested.
  if (nargout >= 3)
    % return the (fairly pointless) uniform bin widths.
    widths = ones(size(abnew)) .* w;
  end
end

