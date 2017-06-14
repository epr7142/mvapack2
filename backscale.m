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
## @anchor{backscale}
## @deftypefn {Function File} {@var{B} =} backscale (@var{A}, @var{center}, @var{scale})
## @deftypefnx {Function File} {@var{B} =} backscale (@var{mdl})
## @deftypefnx {Function File} {[@var{B}, @var{m}] =} backscale (@var{mdl})
## @deftypefnx {Function File} {[@var{B}, @var{m}, @var{s}] =} backscale (@var{mdl})
## Undo a scaling operation performed while building a model of the matrix
## @var{A}. The columns (variables) of @var{A} must match the lengths of
## the vectors @var{center} and @var{scale}. The inputs may either be
## a complete set of data to perform backscaling or a PCA, PLS or OPLS
## model. In the latter case, the extracted centering and scaling vectors
## used during backscaling may be optionally returned as @var{m} and @var{s}
## as well.
## @end deftypefn

function [B, m, s] = backscale (A, center, scale)
  % check the number of arguments.
  if (nargin == 1 && isstruct(A))
    % it is assumed that the user has passed a model and wishes to backscale
    % the loadings matrix into the original high-dimensional space.
    mdl = A;

    % check for the model type field.
    if (!isfield(mdl, 'type'))
      % model doesn't have the type field.
      error('backscale: model is missing a type field, or it is invalid');
    end

    % check the model type.
    if (strcmp(mdl.type, 'pca') || strcmp(mdl.type, 'pls'))
      % extract the loadings.
      A = mdl.P';
    elseif (strcmp(mdl.type, 'opls'))
      % extract the predictive loadings.
      A = mdl.Pp';
    else
      % throw an exception.
      error('backscale: unknown or un-handleable model type');
    end

    % extract the scaling results from the X-matrix scaling function.
    center = mdl.mean.X;
    scale = mdl.scale.X;

    % check if the model is of the multiblock variety.
    if (isfield(mdl, 'weight') && isfield(mdl.weight, 'X'))
      % it is multiblock. add in the superscaling values.
      scale = scale .* mdl.weight.X;
    end
  elseif (nargin != 3)
    % output an error message.
    print_usage();
  end

  % check if the input matrix needs transposition.
  if (columns(A) != length(center))
    % does it need to be flipped?
    if (rows(A) == length(center))
      % flip the matrix.
      A = A';
    else
      % even flipping won't help.
      error('backscale: input data size does not match center/scale');
    end
  end

  % reshape the centering vector into a row vector.
  center = reshape(center, 1, length(center));

  % backscale the matrix A based on the centering matrix and scaling vector
  % to 'reverse' the effects of the original operations.
  B = ones(rows(A), 1) * center + A * diag(scale);

  % determine whether a second output was requested.
  if (nargout >= 2)
    % yes. return the centering vector.
    m = center;
  end

  % determine whether a third output was requested.
  if (nargout >= 3)
    % yes. return the scaling vector.
    s = scale;
  end
end

