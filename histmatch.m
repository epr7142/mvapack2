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
## @anchor{histmatch}
## @deftypefn {Function File} {@var{Xn} =} histmatch (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{s}] =} histmatch (@var{X})
## Normalize the observations of a data matrix using the Histogram
## Matching method described in:
##
## @quotation
## R. J. O. Torgrip et. al. `A note on normalization of biofluid 1H-NMR data'.
## Metabolomics, 2008. 4(2008):114-121.
## @end quotation
## @end deftypefn

function [Xnew, s] = histmatch (X)
  % check for proper arguments.
  if (nargin != 1 || !any(nargout == [1 : 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the input matrix is complex.
  if (iscomplex(X))
    % only use the real portion.
    X = real(X);
  end

  % build a log-transformed vector of data matrix intensities.
  Xall = reshape(X, rows(X) * columns(X), 1);
  Zall = log2(abs(Xall) + 1);

  % linearly map the transformed intensities.
  Zmin = min(Zall);
  Zmax = max(Zall);
  Zmap = linspace(Zmin, Zmax, 50);

  % build a histogram of the median spectrum.
  Xt = median(X);
  Zt = log2(abs(Xt) + 1);
  Ht = hist(Zt, Zmap);

  % initialize the output matrix.
  Xnew = zeros(size(X));
  stmp = ones(rows(X), 1);

  % loop through the rows of the data matrix.
  for i = 1 : rows(X)
    % extract the data matrix row.
    Xs = X(i,:);

    % normalize the data matrix row based on its histogram.
    [a, Fa, info] = fminbnd(@(x) histmatch_func(x, Ht, Xs, Zmap), 0.05, 20);

    % store the normalized data matrix row.
    Xnew(i,:) = X(i,:) .* a;
    stmp(i) = 1 ./ a;
  end

  % check if the normalization values were requested.
  if (nargout == 2)
    % return the requested values.
    s = stmp;
  end
end

