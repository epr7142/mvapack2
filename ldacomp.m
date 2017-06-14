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
## @anchor{ldacomp}
## @deftypefn {Function File} {[@var{P}, @var{D}, @var{U}] =} ldacomp (@var{X}, @var{Y})
## Extracts a multiclass LDA decomposition from a data and response matrix.
## This function is not to be used directly; it is a subroutine of @ref{lda}.
## @end deftypefn

function [P, D, U] = ldacomp (X, Y)
  % get the data matrix sizes.
  [N, K] = size(X);
  M = columns(Y);

  % initialize the covariance matrix partitions.
  Sb = zeros(K, K);
  Sw = zeros(K, K);

  % initialize the means matrix.
  U = [];

  % loop through the classes to calculate the covariance contributions.
  for m = 1 : M
    % store the class indices for extracting observations by class.
    idx = classidx(Y, m);
    n = length(idx);

    % store the data values for the current class.
    Xm = X(idx,:);
    Xc = center(Xm);

    % store the mean observation value for the current class.
    mu = mean(Xm);
    U = [U; mu];

    % calculate the current class' contribution to between-class variation.
    Sb += n .* mu' * mu;

    % calculate the current class' contribution to within-class variation.
    Sw += Xc' * Xc;
  end

  % calculate the eigendecomposition of the covariance matrix quotient.
  [P, D] = eigsort(Sb, Sw);
end

