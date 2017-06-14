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
## @anchor{binoptim}
## @deftypefn {Function File} {@var{xnew} =} binoptim (@var{X}, @var{ab}, @var{w}, @var{slack})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}] =} binoptim (@var{X}, @var{ab}, @var{w}, @var{slack})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}] =} binoptim (@var{X}, @var{ab}, @var{w}, @var{slack})
## @deftypefnx {Function File} {[@var{xnew}, @var{abnew}, @var{widths}, @var{indices}] =} binoptim (@var{X}, @var{ab}, @var{w}, @var{slack})
## Optimally bin a one-dimensional spectrum or spectral dataset in @var{X}
## such that final bins have a width no greater than @var{w}. The optionally
## returnable values in @var{abnew} correspond to the new bin centers in
## abscissa units. The final optional return value (@var{widths}) provides
## the widths of all the output bins.
##
## This code is based on the description of the Optimized Bucketing
## Algorithm (OBA) presented in:
##
## @quotation
## Sousa et. al., `Optimized Bucketing of NMR spectra: Three case studies',
## Chemometrics and Intelligent Lab Systems, 2013.
## @end quotation
##
## The @var{w} value is optional and has a default value of @var{0.025}. Also
## optional is the @var{slack} variable, which must range from @var{0} to
## @var{1} and has a default value of @var{0.45}.
## @end deftypefn

function [xnew, abnew, widths, indices] = binoptim (X, ab, w, slack)
  % check for proper arguments.
  if (!any(nargin == [2 : 4]) || !any(nargout == [1 : 4]) || ...
      !ismatrix(X) || !isvector(ab))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('binoptim: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % check if a third argument was supplied.
  if (nargin < 3)
    % no. use the default initial bin width of 0.025 abscissa units.
    w = 0.025;
  end

  % check if a fourth argument was supplied.
  if (nargin < 4)
    % no. use the default slack value of 0.45.
    slack = 0.45;
  end

  % initialize variables using during binning.
  kbin = floor(K / floor(range(ab) / w));
  s = floor(slack * kbin);
  xbar = mean(X);

  % initialize output variables.
  v = [];
  xnew = [];
  abnew = [];
  widths = [];

  % find the jump locations in the abscissa.
  abjump = findjumps(ab) - 1;
  kstops = abjump + 1;

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

    % store the bin boundary and increment the loop counter.
    v = [v; k1];
    k = k2 + 1;
  end

  % finalize the bin boundary vector.
  v = [v; K];

  % loop through the initial bin boundaries.
  for i = 2 : length(v) - 1
    % avoid jumps.
    if (any(kstops == v(i)))
      continue
    end

    % replace the boundary with the local minimum that exists in the slack.
    xmin = min(xbar(max([v(i)-s; v(i-1)+1]) : min([v(i)+s; v(i+1)-1])));
    kmin = find(xbar == xmin)(1);
    v(i) = kmin;
  end

  % initialize the final bin indices.
  indices = [];

  % loop through the optimized bin boundaries.
  for i = 1 : length(v) - 1
    % get the bin start and end indices.
    k1 = v(i);
    k2 = v(i + 1) - 1;

    % avoid jumps and zero-length bins.
    if (k1 >= k2)
      continue
    end

    % build the segments over which to calculate the next bin.
    Xcut = X(:,k1 : k2);
    abcut = ab(k1 : k2);

    % build the next output bin, abscissa center, width and index.
    xnew = [xnew, abs(trapz(abcut, Xcut, 2))];
    abnew = [abnew; median(abcut)];
    widths = [widths; abs(range(abcut))];
    indices = [indices; k1, k2];
  end
end

