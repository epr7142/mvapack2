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
## @anchor{rqinfo}
## @deftypefn {Function File} {} rqinfo (@var{mdl})
## Prints @math{R^2}/@math{Q^2} values from a PCA, PLS, OPLS or LDA model.
## @end deftypefn

function rqinfo (mdl)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % get the cumulative and component statistics.
  rq_comp = rq(mdl);
  rq_cum = rq(mdl);

  % print a header.
  printf('\n#\tR2\tR2(cum)\tQ2\tQ2(cum)\n');

  % loop for every componennt.
  for a = 1 : mdl.A
    % print the current component statistics.
    printf('%d:\t%.3f\t%.3f\t%.3f\t%.3f\n', a, ...
           rq_comp(a,1), rq_comp(a,2), rq_cum(a,1), rq_cum(a,2));
  end

  % print a footer.
  printf('\n');
end

