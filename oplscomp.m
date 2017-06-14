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
## @anchor{oplscomp}
## @deftypefn {Function File} {[@var{w}, @var{t}, @var{p}, @var{u}, @var{c}, @var{Wo}, @var{To}, @var{Po}, @var{iter}] =} oplscomp (@var{X}, @var{Y}, @var{V}, @var{aout})
## Extracts a single OPLS component from a data and response matrix, as well
## as any accompanying significant orthogonal components. The returned
## @var{w}, @var{t}, @var{p}, @var{u} and @var{c} are vectors corresponding
## to the rank-one approximation of @var{X} and @var{Y} of the OPLS component.
## The OPLS components may be returned in batches of more than one at a time.
## The final input argument @var{aout} may be set to nonzero to specify a
## desired orthogonal component count.
## @end deftypefn

function [w, t, p, u, c, Wo, To, Po, iter] = oplscomp (X, Y, V, aout)
  % initialize the inputs.
  reiterate = true;
  [N, K] = size(X);
  M = columns(Y);
  Ao = 0;
  E = X;

  % calculate the component limit.
  Alim = min([N / 2; K / 2]);

  % initialize the outputs.
  iters = [];
  Wo = [];
  To = [];
  Po = [];

  % loop until the finalization criterion is met.
  while (reiterate == true)
    % calculate a new component.
    [w, t, p, u, c, iter] = plscomp (E, Y);
    iters = [iters; iter];

    % calculate the orthogonal loadings, weights.
    po = p;
    for m = 1 : M
      v = V(:,m);
      po = po - ((v' * po) ./ (v' * v)) * v;
    end
    worth = po;

    % normalize the loadings, weights.
    wo = po ./ sqrt(po' * po);
    to = E * wo;
    po = (E' * to) ./ (to' * to);

    % see if the orthogonal component is significant.
    ev = sqrt(worth' * worth) / sqrt(p' * p);
    if ((aout && Ao < aout) || ...
        (!aout && (!Ao || (ev >= 0.01 && Ao < Alim))))
      % it is significant. store the component.
      Wo = [Wo, wo];
      To = [To, to];
      Po = [Po, po];

      % subtract the orthogonal contribution.
      E = E - to * po';
      Ao++;
    else
      % stop iterating.
      reiterate = false;
    end
  end

  % make sure we found significant orthogonal variation.
  if (columns(Wo) < 1)
    % we didn't. should this be acceptable and handled differently?
    error('oplscomp: failed to extract any significant orthogonal component');
  end

  % return the total number of iterations.
  iter = sum(iters);
end

