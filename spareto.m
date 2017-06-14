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
## @anchor{spareto}
## @tex
## Pareto scaling is a middle-ground between applying no scaling and applying
## UV scaling. Large variations are still weighted more strongly than smaller
## variations, but they will be less likely to dominate the resulting models
## as would occur in unscaled data.
## @end tex
## @deftypefn {Function File} {@var{Xs} =} spareto (@var{X})
## @deftypefnx {Function File} {[@var{Xs}, @var{mu}] =} spareto (@var{X})
## @deftypefnx {Function File} {[@var{Xs}, @var{mu}, @var{s}] =} spareto (@var{X})
## @deftypefnx {Function File} {@dots{} =} spareto (@var{X}, @var{w})
## Scale the variables of a data matrix by the square root of their
## sample standard deviation. Centering is also performed in the
## process. An optional weighting vector @var{w} may be passed during the
## scaling. The variables used to center and scale @var{X} may be
## optionally returned.
## @tex
##
## The resulting scaled elements of $\tilde{X}$ are calculated as follows:
## $$ \tilde{x}_{ik} = { x_{ik} - \bar{x}_k \over w_k \sqrt{s_k} } $$
## @end tex
## @end deftypefn

function [Xnew, mu, s] = spareto (X, w)
  % check for proper arguments.
  if (!any(nargin == [1, 2]) || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a weight vector was passed.
  if (nargin < 2 || !isvector(w) || length(w) != columns(X))
    % no. build a unit weight vector.
    w = ones(1, columns(X));
  else
    % reshape the weight vector.
    w = reshape(w, 1, columns(X));
  end

  % return the scaled output matrix.
  Xnew = center(X) * diag(1 ./ (w .* sqrt(std(X))));

  % see if the center was requested.
  if (nargout >= 2)
    % yes. return the center.
    mu = mean(X);
  end

  % see if the scale was requested.
  if (nargout >= 3)
    % yes. return the scale.
    s = sqrt(std(X));
  end
end

