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
## @anchor{integrals}
## @deftypefn {Function File} {@var{Imax} =} integrals (@var{X}, @var{ab}, @var{roi})
## @deftypefnx {Function File} {[@var{I}, @var{Iab}] =} integrals (@var{X}, @var{ab}, @var{roi})
## Calculates integrals of a spectral dataset over the specified regions of
## interest in @var{roi}. If a single output is requested (@var{Imax}), it will
## hold the final values of the integals. If two outputs are requested
## (@var{I}, @var{Iab}), they will hold the integration curves.
## @end deftypefn

function [Ia, Ib] = integrals (X, ab, roi)
  % check for proper arguments.
  if (nargin != 3 || !any(nargout == [1 : 2]) || ...
      !ismatrix(X) || !isvector(ab) || ...
      !ismatrix(roi) || columns(roi) != 2)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % only use the real portion of the data matrix.
  if (iscomplex(X))
    % keep it real.
    X = real(X);
  end

  % check if the input matrix is actually a vector (single observation case).
  if (isvector(X))
    % reshape into a row vector to obey the matrix rules.
    X = reshape(X, 1, length(X));
  end

  % get the size of the input matrix.
  [N, K] = size(X);

  % check that the length of the abscissa matches the number of variables.
  if (length(ab) != K)
    % no match. throw an exception.
    error('integrals: variable count and abscissa length do not match');
  end

  % reshape the abscissa into a column vector.
  ab = reshape(ab, K, 1);

  % build the boundary matrix from the roi matrix.
  v = sort(findnearest(ab, roi), 2);

  % initialize output variables.
  Iv = [];
  Iab = [];
  Imax = [];

  % loop through the bin boundaries.
  for i = 1 : rows(v)
    % get the integral start and end indices.
    k1 = v(i, 1);
    k2 = v(i, 2) - 1;

    % build the segments over which to integrate.
    Xi = X(:, k1 : k2);
    abi = sort(ab(k1 : k2));

    % build the next output bin, abscissa center, width and index.
    Iv = [Iv, cumtrapz(abi', Xi, 2)];
    Iab = [Iab; abi];
    Imax = [Imax, trapz(abi', Xi, 2)];
  end

  % check if we should columnize the output vector.
  if (isvector(X))
    % yes.
    Iv = Iv';
    Imax = Imax';
  end

  % check the output argument count.
  if (nargout == 1)
    % return just the maxima.
    Ia = Imax;
  elseif (nargout == 2)
    % return the curves.
    Ia = Iv;
    Ib = Iab;
  end
end

