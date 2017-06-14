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
## @anchor{pqnorm}
## @deftypefn {Function File} {@var{Xn} =} pqnorm (@var{X})
## @deftypefnx {Function File} {[@var{Xn}, @var{s}] =} pqnorm (@var{X})
## Normalize the observations of a data matrix using the Probabilistic
## Quotient Normalization method as described in:
##
## @quotation
## F. Dieterle et. al. `Probabilistic Quotient Normalization as Robust
## Method to Account for Dilution of Complex Biological Mixtures.
## Application in 1H NMR Metabonomics.' Analytical Chemistry, 2006.
## @end quotation
## @end deftypefn

function [Xnew, s] = pqnorm (X)
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

  % constant-sum normalize the data to unit integral.
  [Xcs, Scs] = csnorm(X);

  % calculate the median of the normalized data.
  Xref = median(Xcs);

  % initialize the output data.
  Xnew = zeros(size(X));
  Spq = zeros(rows(X), 1);

  % loop through the observations.
  for i = 1 : rows(X)
    % calculate the ratio of the current observation to the reference.
    q = Xcs(i,:) ./ Xref;

    % remove infinities from the ratio.
    q(isinf(q)) = [];
    q(isnan(q)) = [];

    % calculate the probabilistic quotient.
    qmed = median(q);

    % store the normalized observation.
    Xnew(i,:) = Xcs(i,:) ./ qmed;
    Spq(i) = qmed;
  end

  % check if the second argument was requested.
  if (nargout >= 2)
    % return the requested normalization factors.
    s = Scs .* Spq;
  end
end

