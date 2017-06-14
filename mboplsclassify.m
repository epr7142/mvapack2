## Copyright (C) 2015 University of Nebraska Board of Regents.
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
## @anchor{mboplsclassify}
## @deftypefn {Function File} {@var{Y} =} mboplsclassify (@var{mdl}, @var{X})
## @deftypefnx {Function File} {[@var{Y}, @var{T}] =} mboplsclassify (@var{mdl}, @var{X})
## Predicts responses @var{Y} from one or more observations @var{X} based on
## the MBOPLS model provided in @var{mdl}. The observations in @var{X} are
## transformed by the regression coefficients (@var{B}) and classified
## based on sum of squares to the model classes.
##
## @strong{NOTE:} this function is not meant to be used directly. If you
## want to use a PLS model to classify new observations, use @ref{classify}.
## @end deftypefn

function [Y, T] = mboplsclassify (mdl, X)
  % check the type of arguments.
  if (nargin != 2 || !isstruct(mdl) || !iscell(X))
    % invalid arguments. throw an exception.
    print_usage();
  end

  % ensure the input cell array contains the correct number of blocks.
  if (isempty(X) || length(X) != mdl.B)
    % this won't fly. throw an exception.
    error('mboplsclassify: input data array block count mismatch');
  end

  % ensure the input cell array contains only matrices.
  if (!all(cellfun(@(Xb) ismatrix(Xb), X)))
    % this won't fly either. throw an exception.
    error('mboplsclassify: input data array may contain only matrices');
  end

  % get the number of blocks in the dataset.
  B = length(X);

  % get the number of observations in the dataset.
  N = unique(cellfun(@(Xb) rows(Xb), X));

  % ensure only one row count exists in the data.
  if (!isscalar(N))
    % whoops, this really won't fly. throw an exception.
    error('mboplsclassify: blocks must have the same observation count');
  end

  % get the number of variables in each block.
  K = cellfun(@(Xb) columns(Xb), X);
  Kmod = reshape(cellfun(@(blk) blk.K, mdl.blocks), size(K));

  % ensure the variable counts match that of the model.
  if (!all(K == Kmod))
    % yeah... this won't fly. throw an exception.
    error('mboplsclassify: block variable counts must match the model');
  end

  % build the supermatrix and scale it by the superscaling vector.
  Xsup = cell2mat(reshape(X, 1, B)) * diag(1 ./ mdl.superscale.X);

  % execute the supermatrix classification routine.
  [Y, T] = oplsclassify(mdl, Xsup);
end

