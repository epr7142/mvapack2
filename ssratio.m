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
## @anchor{ssratio}
## @deftypefn {Function File} {@var{r} =} ssratio (@var{a}, @var{b})
## Calculates the ratio of the row sum of squares of @var{a} and
## @var{b}.
## @end deftypefn

function r = ssratio (a, b)
  % easy peasy.
  r = sumsq(vec(a)) / sumsq(vec(b));
end

