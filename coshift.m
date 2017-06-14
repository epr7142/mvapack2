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
## @anchor{coshift}
## @deftypefn {Function File} {@var{Xc} =} coshift (@var{X})
## @deftypefnx {Function File} {[@var{Xc}, @var{lags}] =} coshift (@var{X})
## Performs whole-spectrum correlation-optimized shifting to align peaks in NMR
## spectra to a common target, the average of the dataset. The @var{Xc} output
## is the set of aligned spectra, and @var{X} is the input (real) data matrix
## to be shifted. Optionally, the number of points each spectrum was shifted
## may be returned in @var{lags}.
##
## If individual peaks are misaligned within spectra such that this function
## cannot correct them, use @ref{icoshift}. However, this function will not
## produce warping artifacts.
## @end deftypefn

function [Xc, lags] = coshift (X)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 2]) || !ismatrix(X))
    % print the usage statement.
    print_usage();
  end

  % ensure the input data is real.
  if (iscomplex(X))
    % use only real data.
    X = real(X);
  end

  % get the dimensions of the input data matrix.
  [N, K] = size(X);

  % calculate the mean pseudospectrum of the input data matrix.
  Xhat = mean(X);
  Xc = X;

  % see if the lags were requested.
  if (nargout >= 2)
    % initialize the lag vector.
    lags = zeros(N, 1);
  end

  % loop through each row of the input data matrix.
  for i = 1 : N
    % find the shift value that produces maximal cross-correlation.
    [R, lag] = xcorr(Xhat, X(i,:));
    maxlag = lag(find(R == max(R)));

    % circularly shift the spectrum to the optimal lag value.
    Xc(i,:) = shift(X(i,:), maxlag);

    % determine the direction of the circular shift.
    if (lag > 1)
      % right shift. extend the leftmost value to the far left.
      Xc(i, 1 : lag) = ones(1, lag) .* Xc(i, lag + 1);
    elseif (lag < -1)
      % left shift. extend the rightmost value to the far right.
      Xc(i, end + lag + 1 : end) = ones(1, abs(lag)) .* Xc(i, lag - 1);
    end

    % see if the lag should be returned.
    if (nargout >= 2)
      % store the optimal lag value.
      lags(i) = maxlag;
    end
  end
end

