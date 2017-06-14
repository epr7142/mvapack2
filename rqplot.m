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
## @anchor{rqplot}
## @deftypefn {Function File} {} rqplot (@var{mdl})
## @deftypefnx {Function File} {@var{h} =} rqplot (@var{mdl})
## Builds a bar plot of @math{R^2}/@math{Q^2} values from a PCA, PLS,
## OPLS or LDA model.
## @end deftypefn

function h = rqplot (mdl)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % read the statistics from the model and append a row of zeros to ensure
  % the plotting routine can figure out we have two bars per component.
  v = [rq(mdl); 0, 0];

  % initialize the figure window.
  figure();
  htmp = bar(v);
  title('Model Summary: R^2 vs Q^2');
  xlabel('Components');
  ylabel('R^2, Q^2');

  % set the two bar colors.
  set(htmp(1), 'facecolor', [0, 0, 1]);
  set(htmp(2), 'facecolor', [0, 1, 0]);

  % set the x-axis range.
  set(gca(), 'xlim', [0.5, rows(v) - 0.5]);

  % see if the user requested the figure handle.
  if (nargout == 1)
    % return the handle.
    h = htmp;
  end
end

