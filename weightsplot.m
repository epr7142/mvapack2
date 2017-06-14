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
## @anchor{weightsplot}
## @deftypefn {Function File} {} weightsplot (@var{mdl})
## @deftypefnx {Function File} {} weightsplot (@var{mdl}, @var{coloring})
## @deftypefnx {Function File} {@var{h} =} weightsplot (@var{mdl})
## @deftypefnx {Function File} {@var{h} =} weightsplot (@var{mdl}, @var{coloring})
## Builds a loadings plot from PCA, PLS or OPLS modeled data. An optional
## second argument @var{coloring} may be specified to color the points.
## See @ref{nocolors}, @ref{varcolors}.
## @end deftypefn

function h = weightsplot (mdl, coloring)
  % check for a minimum number of arguments.
  if (nargin < 1 || !isstruct(mdl))
    % print the usage statement.
    print_usage();
  end

  % check that the model has enough dimensions.
  if (mdl.A < 2)
    % nope. throw an exception.
    error('weightsplot: at least two model dimensions required');
  end

  % extract the weights.
  W = weights(mdl, 2);
  w1 = W(:,1);
  w2 = W(:,2);

  % see if the coloring was specified.
  if (nargin >= 2 && coloring == @varcolors)
    % build the coloring matrix.
    colors = coloring(mdl.X);
  else
    % make a black matrix.
    colors = [0, 0, 0];
  end

  % initialize the figure.
  figure();
  hold on;
  title('Weights plot: w_1 vs w_2');
  xlabel('w_1');
  ylabel('w_2');
  htmp = scatter(w1, w2, [], colors, '^');

  % plot the main weights ellipse.
  E = ellipse(W, false);
  plot(E.xy(:,1), E.xy(:,2), 'color', [0, 0, 0], 'linewidth', 2, '-');

  % see if the model is a pls model.
  if (strcmp(mdl.type, 'pls'))
    % plot the y weights too.
    c1 = mdl.C(:,1);
    c2 = mdl.C(:,2);
    plot(c1, c2, 'color', [0, 0, 0], '^');
  end

  % finish the plot.
  hold off;

  % see if the plot handle was requested.
  if (nargout == 1)
    % yes. return the plot handle.
    h = htmp;
  end
end

