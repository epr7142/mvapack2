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
## @anchor{icoshift}
## @deftypefn {Function File} {@var{Xc} =} icoshift (@var{X}, @var{ab})
## @deftypefnx {Function File} {@var{Xc} =} icoshift (@var{X}, @var{ab}, @var{seg})
## @deftypefnx {Function File} {[@var{Xc}, @var{lags}] =} icoshift (@var{X}, @var{ab})
## @deftypefnx {Function File} {[@var{Xc}, @var{lags}] =} icoshift (@var{X}, @var{ab}, @var{seg})
## @deftypefnx {Function File} {[@var{Xc}, @var{lags}] =} icoshift (@var{X}, @var{ab}, @var{seg}, @var{cofirst})
## Performs interval correlation-optimized shifting to align peaks in NMR
## spectra to a common target, by default the average of the dataset.
##
## The @var{Xc} output is the set of aligned spectra, and @var{X} is the input
## (real) data matrix to be shifted. The second argument @var{ab} is the
## abscissa of the spectral data.
##
## The optional argument @var{seg} may be used to either define a number of
## segments to use (scalar) or define one or more segments based on abscissa
## values (matrix). Each manually defined segment should be a two-element row
## in @var{seg}. The final optional argument @var{cofirst} allows the user to
## perform an initial global alignment (@ref{coshift}) prior to segmented
## alignment. The default value of @var{cofirst} is @code{false}.
##
## This code is based on the algorithm documented in:
##
## @quotation
## F. Savorani, G. Tomasi, S. B. Engelsen, `icoshift: A versatile tool for
## the rapid alignment of 1D NMR spectra.', J. Magn. Res. 2010: 190-202.
## @end quotation
## @end deftypefn

function [Xc, lags] = icoshift (X, ab, seg, cofirst = false)
  % ensure the minimum number of arguments is met.
  if (!any(nargin == [2 : 4]) || !any(nargout == [1 : 2]))
    % nope. throw a hissy fit.
    print_usage();
  end

  % check the input data matrix.
  if (!ismatrix(X) || isvector(X))
    % invalid type. throw an exception.
    error('icoshift: input data matrix must have at least two observations');
  end

  % check the input abscissa vector.
  if (!isvector(ab) || !isreal(ab))
    % invalid type. throw an exception.
    error('icoshift: abscissa argument must be a real vector');
  end

  % ensure the input data is real.
  if (iscomplex(X))
    % use only real data.
    X = real(X);
  end

  % get the dimensions of the input data matrix.
  [N, K] = size(X);

  % check the length of the abscissa.
  if (length(ab) != K)
    % yikes! throw an exception.
    error(['icoshift: abscissa vector length and data matrix column ', ...
           'count must match']);
  end

  % see if an initial coshift is requested.
  if (cofirst == true)
    % yes. use global correlation-optimized shifting to initialize the input
    % data matrix for further segmented shifting. this minimizes the
    % amount of required segmented shifting.
    Xc = coshift(X);
  else
    % no. use the input matrix as-is.
    Xc = X;
  end

  % calculate the mean pseudospectrum of the input data matrix.
  Xhat = mean(Xc);

  % see if the user neglected the specification of how many segments.
  if (nargin < 3 || isempty(seg))
    % default to fifty segments, found by optimized binning.
    seg = 50;
  end

  % check the type of the segment argument.
  if (isscalar(seg))
    % initialize the segment boundaries using optimized binning.
    [B, abB, wB, idxB] = binoptim(Xc, ab, range(ab) / seg, 0.25);

    % the number of segments is the number of bins.
    M = columns(B);
  elseif (ismatrix(seg) && columns(seg) == 2)
    % build the index matrix based on user-defined segment boundaries.
    idxB = sort(findnearest(ab, seg), 2);

    % the number of segments is the the number of segment matrix rows.
    M = rows(seg);
  else
    % of all the possible argument types, the user has managed to fail.
    error('icoshift: invalid segment argument type');
  end

  % calculate the segment point counts.
  idxB = [idxB, range(idxB, 2) + 1];

  % initialize the lag matrix.
  maxlag = zeros(N, M);

  % loop through each row of the input data matrix.
  for i = 1 : N
    % loop through each segment of the data matrix.
    for j = 1 : M
      % build variables that hold the current segment boundaries.
      p = idxB(j,1);
      q = idxB(j,2);

      % extract the current segments from the data matrix.
      ab_ft = ab(p : q);
      X_ft = Xc(i, p : q);
      Xhat_ft = Xhat(p : q);

      % find the shift value that produces maximal cross-correlation.
      [R, lag] = xcorr(Xhat_ft, X_ft);
      maxlag(i,j) = lag(find(R == max(R)));
      L = maxlag(i,j);

      % limit the shift to a quarter of the segment width.
      if (abs(L) > idxB(j,3) / 4)
        % too far. don't do a shift.
        L = 0;
      end

      % circularly shift the segment to the optimal lag value.
      Xshift = shift(X_ft, L);

      % determine the direction of the circular shift.
      if (L > 0)
        % right shift.
        if (j > 1)
          % interpolate the leftmost points from the previous segment.
          x = ab([p - 2, p - 1, p + L, p + L + 1])';
          y = Xc(i, [p - 2, p - 1, p + L, p + L + 1]);
          xi = ab(p : p + L - 1)';
          yi = interp1(x, y, xi, 'pchip');

          % insert the interpolated points into the shifted segment.
          Xshift(1 : L) = yi;
        else
          % extend the leftmost points.
          Xshift(1 : L) = ones(1, L) .* X_ft(1);
        end
      elseif (L < 0)
        % left shift.
        if (j < M)
          % interpolate the rightmost points from the next segment.
          x = ab([q + L - 1, q + L, q + 1, q + 2])';
          y = Xc(i, [q + L - 1, q + L, q + 1, q + 2]);
          xi = ab(q + L + 1 : q)';
          yi = interp1(x, y, xi, 'pchip');

          % insert the interpolated points into the shifted segment.
          Xshift(end + L + 1 : end) = yi;
        else
          % extend the rightmost points.
          Xshift(end + L + 1 : end) = ones(1, abs(L)) .* X_ft(end);
        end
      end

      % store the processed segment back into the data matrix.
      Xc(i, p : q) = Xshift;
    end
  end

  % see if the lags should be returned.
  if (nargout >= 2)
    % return the calculated lag matrix.
    lags = maxlag;
  end
end

