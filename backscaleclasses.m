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
## @anchor{backscaleclasses}
## @deftypefn {Function File} {@var{Y} =} backscaleclasses (@var{mdl})
## Recover the original discriminant analysis class matrix @var{Y}
## @pxref{classes} from a scaled, centered version using @ref{backscale}.
## @end deftypefn

function Y = backscaleclasses (mdl)
  % check for proper arguments.
  if (nargin != 1 || nargout != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check for the required backscaling information.
  if (!isfield(mdl, 'Y') || ...
      !isfield(mdl, 'mean') || !isfield(mdl.mean, 'Y') || ...
      !isfield(mdl, 'scale') || !isfield(mdl.scale, 'Y'))
    % throw an exception.
    error('backscaleclasses: insufficient information to rebuild Y');
  end

  % check that the class matrix has at least two columns.
  if (columns(mdl.Y) < 2)
    % throw an exception. class matrices are two or more columns.
    error('backscaleclasses: Y doesnt look like a class matrix');
  end

  % rebuild the classes matrix.
  Y = round(abs(backscale(mdl.Y, mdl.mean.Y, mdl.scale.Y)));
end

