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
## @anchor{plscomp}
## @deftypefn {Function File} {[@var{w}, @var{t}, @var{p}, @var{u}, @var{c}, @var{iter}] =} plscomp (@var{X}, @var{Y})
## Extracts a single PLS component from a data and response matrix.
## @end deftypefn

function [w, t, p, u, c, iter] = plscomp (X, Y)
  % initialize the x and y scores.
  t = rand(rows(X), 1);
  u = Y(:,1);

  % initialize the tolerance.
  tol = 1.0e-9;

  % initialize.
  told = t;
  iter = 1;

  % loop until convergence.
  while (iter == 1 || norm(told - t) / norm(t) > tol)
    % save the old score vector.
    told = t;

    % calculate the next estimate of the weights, scores and loadings.
    w = (X' * u) ./ (u' * u);
    w = w ./ sqrt(w' * w);
    t = X * w;
    c = (Y' * t) ./ (t' * t);
    u = (Y * c) ./ (c' * c);

    % increment the iteration count.
    iter++;

    % adaptively grow the tolerance to ensure termination.
    if (mod(iter, 100) == 0)
      tol *= 10;
    end
  end

  % return the x loading.
  p = (X' * t) ./ (t' * t);
end

