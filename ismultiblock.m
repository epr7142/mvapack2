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
## @anchor{ismultiblock}
## @deftypefn {Function File} {@var{tf} =} ismultiblock (@var{mdl})
## Returns whether a model @var{mdl} is of the `multiblock' variety.
## @end deftypefn

function tf = ismultiblock (mdl)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % determine if any jumps were found.
  if (isfield(mdl, 'B') && isfield(mdl, 'blocks') && ...
      isscalar(mdl.B) && mdl.B == length(mdl.blocks))
    % yes! multiblock.
    tf = true;
  else
    % no. uniform.
    tf = false;
  end
end

