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
## @anchor{loadingsplot}
## @deftypefn {Function File} {} loadingsplot (@var{mdl})
## @deftypefnx {Function File} {} loadingsplot (@var{mdl}, @var{coloring})
## @deftypefnx {Function File} {} loadingsplot (@var{mdl}, @var{coloring}, @var{numbers})
## @deftypefnx {Function File} {@var{h} =} loadingsplot (@var{mdl})
## @deftypefnx {Function File} {@var{h} =} loadingsplot (@var{mdl}, @var{coloring})
## @deftypefnx {Function File} {@var{h} =} loadingsplot (@var{mdl}, @var{coloring}, @var{numbers})
## Builds a loadings plot from PCA, PLS or OPLS modeled data. An optional
## second argument @var{coloring} may be set to an appropriate function
## (See @ref{nocolors}, @ref{varcolors}). An optional third argument
## (@var{numbers}) may be set to true (default is false) to change the
## plotted points into variable numbers.
## @end deftypefn

function h = loadingsplot (mdl, coloring, numbers)
  % check for proper arguments.
  if (!any(nargin == [1 : 3]) || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the model contains enough components.
  if (mdl.A < 2)
    % nope. throw an exception.
    error('loadingsplot: at least two model dimensions required');
  end

  % extract the loadings from the model.
  P = loadings(mdl, 2);
  p1 = P(:,1);
  p2 = P(:,2);

  % check if a custom coloring method was selected.
  if (nargin < 2 || isempty(coloring))
    % default to variable coloring.
    coloring = @varcolors;
  else
    % make sure the coloring method is valid.
    if (!is_function_handle(coloring))
      % not valid. throw an exception.
      error('loadingspot: coloring argument must be a valid function handle');
    end
  end

  % build the coloring matrix.
  colors = coloring(mdl.X);

  % check if numbers were requested.
  if (nargin < 3 || isempty(numbers) || !isbool(numbers))
    % default to plotting points.
    numbers = false;
  end

  % initialize the figure.
  figure();
  hold on;
  title("Loadings plot: p_1 vs p_2");
  xlabel("p_1");
  ylabel("p_2");

  % see if numbers or points are desired.
  if (numbers == true)
    % loop through the loadings.
    for k = 1 : rows(P)
      % print the current variable index.
      text(P(k,1), P(k,2), num2str(k), 'color', colors(k,:));
    end
  else
    % scatter plot the score values.
    htmp = scatter(p1, p2, [], colors, '^');
  end

  % line plot the group ellipses.
  E = ellipse(P, false);
  plot(E.xy(:,1), E.xy(:,2), 'color', [0, 0, 0], 'linewidth', 2, '-');

  % release the figure for plotting.
  hold off;

  % see if the figure handle was requested.
  if (nargout == 1)
    % yes. return the handle.
    h = htmp;
  end
end

