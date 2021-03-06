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
## @anchor{mbpcaclassify}
## @deftypefn {Function File} {@var{Y} =} mbpcaclassify (@var{mdl}, @var{X})
## @deftypefnx {Function File} {[@var{Y}, @var{T}] =} mbpcaclassify (@var{mdl}, @var{X})
## Predicts responses @var{Y} from one or more observations @var{X} based on
## the MBPCA model provided in @var{mdl}. The observations in @var{X} are
## transformed into the principal component space and classified based on
## Mahalanobis distances to the model classes.
##
## @strong{NOTE:} this function is not meant to be used directly. If you
## want to use a PCA model to classify new observations, use @ref{classify}.
## @end deftypefn

function [Y, T] = mbpcaclassify (mdl, X)
  % check the type of arguments.
  if (nargin != 2 || !isstruct(mdl) || !iscell(X))
    % invalid arguments. throw an exception.
    print_usage();
  end

  % ensure the input cell array contains the correct number of blocks.
  if (isempty(X) || length(X) != mdl.B)
    % this won't fly. throw an exception.
    error('mbpcaclassify: input data array block count mismatch');
  end

  % ensure the input cell array contains only matrices.
  if (!all(cellfun(@(Xb) ismatrix(Xb), X)))
    % this won't fly either. throw an exception.
    error('mbpcaclassify: input data array may contain only matrices');
  end

  % get the number of blocks in the dataset.
  B = length(X);

  % get the number of observations in the dataset.
  N = unique(cellfun(@(Xb) rows(Xb), X));

  % ensure only one row count exists in the data.
  if (!isscalar(N))
    % whoops, this really won't fly. throw an exception.
    error('mbpcaclassify: blocks must have the same observation count');
  end

  % get the number of variables in each block.
  K = cellfun(@(Xb) columns(Xb), X);
  Kmod = reshape(cellfun(@(blk) blk.K, mdl.blocks), size(K));

  % ensure the variable counts match that of the model.
  if (!all(K == Kmod))
    % yeah... this won't fly. throw an exception.
    error('mbpcaclassify: block variable counts must match the model');
  end

  % build the supermatrix and scale it by the superscaling vector.
  Xsup = cell2mat(reshape(X, 1, B)) * diag(1 ./ mdl.superscale.X);

  % extract the class matrix and scores from the model.
  Ymod = backscaleclasses(mdl);
  Tmod = mdl.T;

  % transform the observations into the PC space.
  Tnew = Xsup * mdl.P;

  % initialize the distance matrix.
  d = zeros(rows(X), mdl.M);

  % loop through the embedded classes in the model.
  for m = 1 : mdl.M
    % extract the scores for the current class.
    idx = classidx(Ymod, m);
    Tcls = Tmod(idx,:);

    % calculate the mean and inverse covariance of the class scores.
    u = mean(Tcls);
    C = inv(cov(Tcls));

    % calculate the mahalanobis distances from the observations
    % to the currently indexed class.
    dv = Tnew - ones(rows(Tnew), 1) * u;
    d(:,m) = diag(dv * C * dv');
  end

  % initialize the output matrix.
  Y = [];
  YM = eye(mdl.M);

  % locate the minimum distance for each observation.
  [dmin, Mmin] = min(d, [], 2);

  % generate a response value for all observations simultaneously.
  Y = YM(Mmin,:);

  % see if a matrix of scores was requested.
  if (nargout >= 2)
    % yes. return the scores.
    T = Tnew;
  end
end

