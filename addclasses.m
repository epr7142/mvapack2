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
## @anchor{addclasses}
## @deftypefn {Function File} {@var{mdlAdd} =} addclasses (@var{mdl}, @var{Y})
## @deftypefnx {Function File} {@var{mdlAdd} =} addclasses (@var{mdl}, @var{Y}, @var{overwrite})
## Adds a supplementary class matrix to a PCA, PLS, @emph{etc} model. This
## function will not work for supervised models that already contain a class
## matrix unless the @var{overwrite} argument is set to @code{true}.
##
## The function is intended to be used like so: @*
## @code{mdl = addclasses(mdl, Y);}
## @end deftypefn

function mdlAdd = addclasses (mdl, Y, overwrite)
  % check for proper arguments.
  if (!any(nargin == [2 : 3]) || nargout != 1 || ...
      !isstruct(mdl) || !ismatrix(Y))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check for a third boolean argument.
  if (nargin < 3 || !isbool(overwrite))
    % default to not overwriting.
    overwrite = false;
  end

  % check if the field already has a class matrix.
  if (isfield(mdl, 'Y') && overwrite == false)
    % it does. user will have to remove it first, then add then new one.
    error('addclasses: model already contains a class matrix');
  end

  % ensure the observation count matches the response count.
  if (mdl.N != rows(Y))
    % whoops. throw an exception.
    error('addclasses: row counts of data and class matrices do not match');
  end

  % return the updated model as a copy of the original.
  mdlAdd = mdl;

  % update the discrimination status.
  mdlAdd.isda = isclasses(Y);

  % add the original and the scaled classes to the model structure.
  mdlAdd.Y0 = Y;
  mdlAdd.M = columns(Y);
  [mdlAdd.Y, mdlAdd.mean.Y, mdlAdd.scale.Y] = mdl.scaling(Y);

  % see if the model is multiblock.
  if (ismultiblock(mdlAdd))
    % add the classes to all blocks as well.
    for b = 1 : mdlAdd.B
      mdlAdd.blocks{b} = addclasses(mdlAdd.blocks{b}, Y, overwrite);
    end
  end
end

