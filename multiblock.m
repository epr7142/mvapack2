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
## @anchor{multiblock}
## @deftypefn {Function File} {@var{X} =} multiblock (@var{X_1}, @dots{}, @var{X_B})
## Builds a multiblock data matrix from individual data matrices, ensuring
## that all dimensions match up.
## @end deftypefn

function X = multiblock (varargin)
  % check for proper arguments.
  if (nargin < 2 || nargout != 1)
    % improper arguments. print the usage statement.
    print_usage();
  end

  % how convenient!
  X = varargin;

  % ensure the input cell array contains only matrices.
  if (!all(cellfun(@(Xb) ismatrix(Xb), X)))
    % this won't fly either. throw an exception.
    error('multiblock: input data array may contain only matrices');
  end

  % ensure only one row count exists in the data.
  if (!isscalar(unique(cellfun(@(Xb) rows(Xb), X))))
    % whoops, this really won't fly. throw an exception.
    error('multiblock: blocks must have the same observation count');
  end
end

