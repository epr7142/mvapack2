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
## @anchor{pcacomp}
## @deftypefn {Function File} {[@var{t}, @var{p}, @var{iter}] =} pcacomp (@var{X})
## @deftypefnx {Function File} {[@var{t}, @var{p}, @var{iter}] =} pcacomp (@var{X}, @var{t0})
## Extracts a single PCA component from a data matrix. The optional
## second argument @var{t0} can be passed to specify an initial value
## for @var{t} prior to iteration.
## @end deftypefn

function [t, p, iter] = pcacomp (X, t0)
  % check if an initial score vector was provided.
  if (nargin == 2 && isvector(t0))
    % yes. use it.
    t = t0;
  else
    % no. use a random vector.
    t = rand(rows(X), 1);
  end

  % initialize.
  told = t;
  iter = 1;

  % loop until convergence.
  while (iter == 1 || norm(told - t) / norm(t) > 1.0e-9)
    % store the previous score.
    told = t;

    % calculate the next estimate of scores and loadings.
    p = (X' * t) ./ (t' * t);
    p = p ./ sqrt(p' * p);
    t = X * p;

    % increment the iteration count.
    iter++;
  end
end

