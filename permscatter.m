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
## @anchor{permscatter}
## @deftypefn {Function File} {} permscatter (@var{S})
## Plots information in @var{S} that has been calculated by @ref{permtest},
## in a scatter plot format.
## @end deftypefn

function permscatter (S)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(S))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % build the legend strings.
  r2str = sprintf('R^2 (p = %.02e)', S.Rsq.p);
  q2str = sprintf('Q^2 (p = %.02e)', S.Qsq.p);
  legstr = {r2str, q2str};

  % build the x-axis label string.
  xstr = sprintf('Response correlation: %d permutations, %d components', ...
                 S.n, S.A);

  % initialize the figure.
  figure();
  hold on;
  title('Permutation test results');
  xlabel(xstr);
  ylabel('R^2 / Q^2');

  % plot the cross-validation lines.
  plot(S.Rsq.xy(:,1), S.Rsq.xy(:,2), 'color', [0, 0, 1]);
  plot(S.Qsq.xy(:,1), S.Qsq.xy(:,2), 'color', [0, 1, 0]);

  % add a plot legend.
  legend(legstr, 'location', 'southeast');

  % plot the cross-validation points.
  scatter([1; S.r], [S.Rsq.orig; S.Rsq.perm], [], [0, 0, 1]);
  scatter([1; S.r], [S.Qsq.orig; S.Qsq.perm], [], [0, 1, 0]);

  % release the figure for plotting.
  hold off;
end

