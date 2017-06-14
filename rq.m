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
## @anchor{rq}
## @deftypefn {Function File} {@var{v} =} rq (@var{mdl})
## Returns the calculated @math{R^2}/@math{Q^2} values from a PCA, PLS
## OPLS or LDA model.
## @end deftypefn

function v = rq (mdl)
  % check for proper arguments.
  if (nargin != 1 || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the requisite fields are available in the model structure.
  if (!isfield(mdl, 'Rsq') || !isfield(mdl, 'Qsq'))
    % invalid model structure: throw an exception.
    error('scores: input model is invalid');
  end

  % see what type of model we're looking at.
  if (strcmp(mdl.type, 'pca'))
    % pca. use the R2X and Q2X.
    v = [cumsum(mdl.Rsq.X.comp), cumsum(mdl.Qsq.comp)];
  elseif (strcmp(mdl.type, 'pls') || ...
          strcmp(mdl.type, 'opls') || ...
          strcmp(mdl.type, 'lda'))
    % pls/opls. use R2Y and Q2Y.
    v = [cumsum(mdl.Rsq.Y.comp), cumsum(mdl.Qsq.comp)];
  else
    % unsupported type. throw an exception.
    error('rq: unsupported model type "%s"', mdl.type);
  end
end

