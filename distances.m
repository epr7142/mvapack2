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
## @anchor{distances}
## @deftypefn {Function File} {@var{D} =} distances (@var{mdl})
## @deftypefnx {Function File} {@var{D} =} distances (@var{mdl}, @var{k})
## @deftypefnx {Function File} {@var{D} =} distances (@var{mdl}, @var{k}, @var{metric})
## @deftypefnx {Function File} {@var{D} =} distances (@var{mdl}, @var{k}, @var{metric}, @var{Y})
## Computes a distance matrix between all class scores of a PCA, PLS or OPLS
## model. PCA models need an accompanying @var{Y}-matrix to define classes.
##
## The default number of components (@var{k}) will be set to the model
## component count, unless specified in the arguments. The default
## metric is the squared Mahalanobis distance.
## @end deftypefn

function D = distances (mdl, k, metric, Y)
  % check for the minimum number of proper arguments: the model.
  if (!any(nargin == [1 : 4]) || !isstruct(mdl))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if the number of components to calculate distances over is provided.
  if (nargin < 2 || isempty(k))
    % not provided. default to the total number of components.
    k = mdl.A;
  end

  % check if the distance metric is provided.
  if (nargin < 3 || isempty(metric))
    % not provided. default to mahalanobis distances.
    metric = @mahalanobis;
  end

  % check if the class matrix is provided.
  if (nargin < 4 || isempty(Y))
    % not provided. check if it exists in the model structure.
    if (!isfield(mdl, 'Y'))
      % nope. throw an exception.
      error('distances: class matrix (Y) required for unsupervised data');
    else
      % no. use the classes embedded in the model.
      Y = backscaleclasses(mdl);
    end
  end

  % extract the desired number of components from the model.
  T = scores(mdl, k);

  % get the number of classes.
  M = columns(Y);

  % initialize the output distance matrix.
  D = zeros(M, M);

  % loop through the first class.
  for m1 = 1 : M
    % get the scores for the first class.
    idx1 = find(Y(:,m1) == 1);
    t1 = T(idx1,:);

    % initialize the diagonal element.
    D(m1,m1) = 0;

    % loop through the second class.
    for m2 = m1 + 1 : M
      % get the scores for the second class.
      idx2 = find(Y(:,m2) == 1);
      t2 = T(idx2,:);

      % calculate the off-diagonal element.
      D(m1,m2) = metric(t1, t2);
      D(m2,m1) = D(m1,m2);
    end
  end
end

