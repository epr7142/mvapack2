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
## @anchor{classidx}
## @deftypefn {Function File} {@var{idx} =} classidx (@var{Y}, @var{m})
## Extracts observation indices of all observations belonging to a given class
## @var{m} from a discriminant analysis class matrix @var{Y}.
## @end deftypefn

function idx = classidx (Y, m)
  % check for proper arguments.
  if (nargin != 2 || !ismatrix(Y) || !isscalar(m) || !isreal(m))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % extract the indices of the observations belonging to the requested class.
  idx = find(Y(:,m) == 1);

  % ensure the number of indices is postitive.
  if (length(idx) <= 0)
    % failure.
    error('classidx: failed to extract class indices: check your Y matrix');
  end
end

