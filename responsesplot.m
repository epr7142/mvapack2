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
## @anchor{responsesplot}
## @deftypefn {Function File} {} responsesplot (@var{mdl})
## @deftypefnx {Function File} {} responsesplot (@var{mdl}, @var{coloring})
## Builds a responses plot from PLS or OPLS modeled data. An optional
## second argument @var{coloring} may be specified to color the points.
## See @ref{nocolors}, @ref{obscolors}.
## @end deftypefn

function h = responsesplot (mdl, coloring)
  % check for a minimum number of arguments.
  if (!any(nargin == [1 : 2]) || !isstruct(mdl))
    % print the usage statement.
    print_usage();
  end

  % extract and backscale the responses.
  [Y, Yfit, Yhat] = responses(mdl);

  % check that the model has enough dimensions.
  if (mdl.M != 1)
    % nope. throw an exception.
    error('responsesplot: function only available for one response variable');
  end

  % see if the coloring was specified.
  if (nargin >= 2 && coloring == @obscolors)
    % build the coloring matrix.
    colors = coloring(mdl.X);
  else
    % make a default color matrix.
    colors = [0, 0, 1];
  end

  % initialize the figure.
  figure();
  hold on;
  title('Responses plot');
  xlabel('y_{in}');
  ylabel('y_{fit} / y_{hat}');
  htmp = scatter(Y, Yhat, [], colors);
  plot(Y, Yfit, 'color', [0, 0, 0], '*');

  % finish the plot.
  hold off;

  % see if the plot handle was requested.
  if (nargout == 1)
    % yes. return the plot handle.
    h = htmp;
  end
end

