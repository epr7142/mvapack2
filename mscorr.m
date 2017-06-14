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
## @anchor{mscorr}
## @deftypefn {Function File} {@var{Xn} =} mscorr (@var{X})
## @deftypefnx {Function File} {@var{Xn} =} mscorr (@var{X}, @var{r})
## Use multiplicative scatter correction to normalize spectra into better
## alignment with each other, row-wise at least. If the optional @var{r}
## reference observation is not provided, the mean of the observations in
## @var{X} will be used.
## @end deftypefn

function [Xnew, b] = mscorr (X, r)
  % check for proper arguments.
  if (!any(nargin == [1 : 2]) || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the input matrix is complex.
  if (iscomplex(X))
    % only use the real portion.
    X = real(X);
  end

  % check if the data is actually a vector.
  if (isvector(X))
    % reshape it into a row vector.
    X = reshape(X, 1, length(X));

    % make sure the reference was also provided.
    if (nargin < 2 || isempty(r))
      % can't have this.
      error(['mscorr: a reference is required for normalizing a ', ...
             'single observation']);
    end
  end

  % check if the reference spectrum was provided.
  if (nargin >= 2 && !isempty(r) && isvector(r))
    % yes. ensure it's properly shaped.
    r = reshape(r, 1, length(r));
  else
    % no. use the default value, the mean.
    r = mean(X);
  end

  % initialize the normalization vector.
  bv = ones(rows(X), 1);

  % initialize the data for normalization.
  Xnew = X;
  rbar = mean(r);
  Xc = center(X, 2);
  rc = (r - rbar)';
  K = inv(rc' * rc) * rc';

  % loop through the rows of the centered data matrix.
  for i = 1 : rows(Xc)
    % normalize and shift the row to best match the target row.
    bv(i) = K * Xc(i,:)';
    Xnew(i,:) = (Xc(i,:)' ./ bv(i) + rbar)';
  end

  % check if a second output was requested.
  if (nargout >= 2)
    % yes. return the normalization factors.
    b = bv;
  end

  % check if the data is actually a vector.
  if (isvector(Xnew))
    % reshape it back into a column vector.
    Xnew = reshape(Xnew, length(Xnew), 1);
  end
end

