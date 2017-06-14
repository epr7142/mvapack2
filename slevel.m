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
## @anchor{slevel}
## @tex
## Level scaling focuses multivariate analysis on fitting the relative
## responses of all signals. Signals with larger relative responses will
## be weighted more heavily than those with smaller relative responses.
## @end tex
## @deftypefn {Function File} {@var{Xs} =} slevel (@var{X})
## @deftypefnx {Function File} {[@var{Xs}, @var{mu}] =} slevel (@var{X})
## @deftypefnx {Function File} {[@var{Xs}, @var{mu}, @var{s}] =} slevel (@var{X})
## @deftypefnx {Function File} {@dots{} =} slevel (@var{X}, @var{w})
## Scale the variables of a data matrix by their mean value,
## the level scaling method. Centering is also performed in the process.
## An optional weighting vector @var{w} may be passed during the scaling.
## The values used to center and scale @var{X} may optionally be returned.
## @tex
##
## The resulting scaled elements of $\tilde{X}$ are calculated as follows:
## $$ \tilde{x}_{ik} = { x_{ik} - \bar{x}_k \over w_k \bar{x}_k } $$
## @end tex
## @end deftypefn

function [Xnew, mu, s] = slevel (X, w)
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

  % calculate the new scaled matrix.
  Xnew = center(X) * diag(1 ./ (w .* mean(X)));

  % see if the center was requested.
  if (nargout >= 2)
    % return the center.
    mu = mean(X);
  end

  % see if the scale was requested.
  if (nargout >= 3)
    % return the scale.
    s = mu;
  end
end

