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
## @anchor{scores}
## @deftypefn {Function File} {@var{T} =} scores (@var{mdl})
## @deftypefnx {Function File} {@var{T} =} scores (@var{mdl}, @var{n})
## Returns the calculated scores values from a PCA, PLS or OPLS model.
## When @var{n} is provided as a second argument, only @var{n} columns
## will be returned in @var{T}. Otherwise, the full matrix of scores
## will be returned.
## @end deftypefn

function T = scores (mdl, n)
  % check for proper arguments.
  if (nargin < 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the requisite fields are available in the model structure.
  if (!isfield(mdl, 'T') || !isfield(mdl, 'A'))
    % invalid model structure: throw an exception.
    error('scores: input model is invalid');
  end

  % check if the component count was supplied.
  if (nargin == 2)
    % extract the desired number of components.
    T = mdl.T(:,1 : min([n; mdl.A]));
  else
    % extract all available components.
    T = mdl.T;
  end
end

