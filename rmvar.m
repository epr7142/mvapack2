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
## @anchor{rmvar}
## @deftypefn {Function File} {[@var{Xrm}, @var{abrm}] =} rmvar (@var{X}, @var{ab}, @var{idx})
## Removes a variable (@var{idx} scalar) or variables (@var{idx} vector)
## from the dataset @var{X} and its corresponding abscissa, @var{ab}.
## Use @ref{findnearest} to locate the indices to supply to @var{idx}.
## @end deftypefn

function [Xrm, abrm] = rmvar (X, ab, idx)
  % check for proper arguments.
  if (nargin != 3 || nargout != 2 || !ismatrix(X) || !isvector(ab) || ...
      !(isvector(idx) || isscalar(idx)))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check for an empty index array.
  if (!isscalar(idx) && !length(idx))
    % output an error.
    error('rmobs: index argument is empty. check your index ordering')
  end

  % initialize the output values.
  Xrm = X;
  abrm = ab;

  % remove the variables from the output values.
  Xrm(:,sort(idx)) = [];
  abrm(sort(idx)) = [];
end

