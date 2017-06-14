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
## @anchor{addlabels}
## @deftypefn {Function File} {@var{mdlAdd} =} addlabels (@var{mdl}, @var{labels})
## Adds a supplementary string cell array (@var{labels}) to a PCA, PLS,
## @emph{etc} model. The labels will be placed on scores plots and written
## to saved scores files. If @var{mdl} already contains a labels array,
## it will be replaced with @var{labels} without warnings.
##
## The function is to be used as follows: @*
## @code{mdl = addlabels(mdl, labels);}
## @end deftypefn

function mdlAdd = addlabels (mdl, labels)
  % check for proper arguments.
  if (nargin != 2 || nargout != 1 || !isstruct(mdl) || !iscellstr(labels))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % return the updated model, with labels added in.
  mdlAdd = mdl;
  mdlAdd.labels = labels;

  % see if the model is multiblock.
  if (ismultiblock(mdlAdd))
    % add the labels to all blocks as well.
    for b = 1 : mdlAdd.B
      mdlAdd.blocks{b} = addlabels(mdlAdd.blocks{b}, labels);
    end
  end
end

