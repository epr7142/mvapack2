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
## @anchor{spassthru}
## @tex
## Performing no scaling prior to multivariate analysis results in fitting
## based on covariance eigenstructure, not correlation eigenstructure. In
## English, large variations will be weighted much more strongly than smaller
## variations in the fitted models.
## @end tex
## @deftypefn {Function File} {@var{Xc} =} spassthru (@var{X})
## @deftypefnx {Function File} {[@var{Xc}, @var{mu}] =} spassthru (@var{X})
## @deftypefnx {Function File} {[@var{Xc}, @var{mu}, @var{s}] =} spassthru (@var{X})
## @deftypefnx {Function File} {@dots{} =} spassthru (@var{X}, @var{w})
## Performs no mean-centering and no scaling. An optional weighting vector
## @var{w} may be passed during the scaling. The variables used to
## center and scale @var{X} may be optionally returned.
## @tex
##
## The resulting scaled elements of $\tilde{X}$ are calculated as follows:
## $$ \tilde{x}_{ik} = { x_{ik} \over w_k } $$
## @end tex
## @end deftypefn

function [Xnew, mu, s] = spassthru (X, w)
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

  % return the centered matrix.
  Xnew = X * diag(1 ./ w);

  % check if the center is requested.
  if (nargout >= 2)
    % return the center.
    mu = zeros(1, columns(X));
  end

  % check if the scale is requested.
  if (nargout >= 3)
    % return the scale.
    s = ones(1, columns(X));
  end
end

