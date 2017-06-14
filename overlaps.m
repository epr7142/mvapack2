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
## @anchor{overlaps}
## @deftypefn {Function File} {@var{D} =} overlaps (@var{mdl}, @var{k})
## @deftypefnx {Function File} {@var{D} =} overlaps (@var{mdl}, @var{k}, @var{Y})
## Computes a p-value matrix between all class scores of a PCA, PLS or OPLS
## model. PCA models need an accompanying @var{Y}-matrix to define classes
## (@xref{addclasses}).
##
## The default number of components (@var{k}) will be set to the model
## component count, unless specified in the arguments.
## @end deftypefn

function p = overlaps (mdl, k, Y)
  % check if the number of expected arguments was passed.
  if (!any(nargin == [1 : 3]) || !isstruct(mdl))
    % print the usage statement.
    print_usage();
  end

  % check if the number of components to calculate distances over is provided.
  if (nargin < 2 || isempty(k))
    % not provided. default to the total number of components.
    k = mdl.A;
  end

  % check if the class matrix is provided.
  if (nargin < 3 || isempty(Y))
    % not provided. check if it should have been provided.
    if (strcmp(mdl.type, 'pca') == 1 && !isfield(mdl, 'Y'))
      % yes. throw an exception.
      error('overlaps: class matrix (Y) required for unsupervised data');
    else
      % no. use the classes embedded in the model.
      Y = backscaleclasses(mdl);
    end
  end

  % get the number of classes.
  M = columns(Y);

  % calculate the distance matrix for the model.
  D = distances(mdl, k, @mahalanobis, Y);

  % initialize the p-value matrix.
  p = ones(M, M);

  % loop through the first class.
  for m1 = 1 : M
    % get the first class indices.
    idx1 = find(Y(:,m1) == 1);
    n1 = length(idx1);

    % set the diagonal element of the output matrix.
    p(m1,m1) = 1;

    % loop through the second class.
    for m2 = m1 + 1 : M
      % get the second class indices.
      idx2 = find(Y(:,m2) == 1);
      n2 = length(idx2);

      % calculate the hotelling T-squared statistic.
      T2 = ((n1 * n2) / (n1 + n2));
      T2 *= ((n1 + n2 - k - 1) / (k * (n1 + n2 - 2)));
      T2 *= D(m1,m2);

      % calculate the p-value from the F distribution.
      p(m1,m2) = 1 - fcdf(T2, k, n1 + n2 - k - 1);
      p(m2,m1) = p(m1,m2);
    end
  end
end

