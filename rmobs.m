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
## @anchor{rmobs}
## @deftypefn {Function File} {@var{Xrm} =} rmobs (@var{X}, @var{idx})
## Removes an observation (@var{idx} scalar) or observations (@var{idx}
## vector) from the dataset @var{X}.
## @end deftypefn

function Xrm = rmobs (X, idx)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1 || !ismatrix(X) || ...
      !(isvector(idx) || isscalar(idx)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check for an empty index array.
  if (!isscalar(idx) && !length(idx))
    % output an error.
    error('rmobs: index argument is empty. check your index ordering')
  end

  % remove the observations from the matrix.
  Xrm = X;
  Xrm(sort(idx),:) = [];
end

