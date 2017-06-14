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
## @anchor{loadings}
## @deftypefn {Function File} {@var{P} =} loadings (@var{mdl})
## @deftypefnx {Function File} {@var{P} =} loadings (@var{mdl}, @var{n})
## Returns the calculated loadings values from a PCA, PLS or OPLS
## model. When @var{n} is provided as a second argument, only @var{n}
## columns will be returned in @var{P}. Otherwise, the full matrix of
## loadings will be returned.
## @end deftypefn

function P = loadings (mdl, n)
  % check for proper arguments.
  if (!any(nargin == [1 : 2]) || nargout != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the requisite fields are available in the model structure.
  if (!isfield(mdl, 'P') || !isfield(mdl, 'A'))
    % invalid model structure: throw an exception.
    error('loadings: input model is invalid');
  end

  % check if the component count was supplied.
  if (nargin == 2)
    % extract the desired number of components.
    P = mdl.P(:, 1 : min([n; mdl.A]));
  else
    % extract all available components.
    P = mdl.P;
  end
end

