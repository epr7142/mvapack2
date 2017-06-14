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
## @anchor{cvindices}
## @deftypefn {Function File} {@var{idx} =} cvindices (@var{N}, @var{Np})
## Creates a vector of group indices from 1 to @var{Np}, where each value
## in the vector is the index of the training set to which the point will
## belong during cross-validation.
## @end deftypefn

function gidx = cvindices (N, Np)
  % generate N random numbers in 1:N, and then
  % restrict them (modulo) to the number of CV groups.
  [p, idx] = sort(rand(N, 1));
  gidx = mod(idx, Np) + 1;
end

