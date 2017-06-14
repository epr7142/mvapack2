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
## @anchor{rqdiff}
## @deftypefn {Function File} {@var{v} =} rqdiff (@var{mdl})
## Returns the calculates @math{R^2}/@math{Q^2} values from a PCA, PLS
## or OPLS model, but gives differential values. In other words, this
## function gives the @math{R^2} and @math{Q^2} values per-component,
## not cumulative.
## @end deftypefn

function v = rqdiff (mdl)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % get the values matrix.
  vcum = rq(mdl);

  % return the differential values matrix.
  if (rows(vcum) == 1)
    % just return the one component of values.
    v = vcum;
  else
    % return the differences.
    v = [vcum(1,:); diff(vcum)];
  end
end

