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
## @anchor{mbpca}
## @tex
## Multiblock Principal Component Analysis (MBPCA) generates a linear
## model of a set of $B$ input data matrices $X = [X_1  X_2 \cdots X_B]$
## such that the overall model components are orthogonal and
## capture a maximal amount of variation present in $X$. The MBPCA
## model is structured as follows:
## $$ X = [X_1@ @  X_2@ @ \cdots@ @ X_B], \forall b \in [1, B] $$
## Where each block ($X_b$) is itself a bilinear model:
## $$ X_b = T_b P_b' + E_b = \sum_{a=1}^A t_{b,a} p_{b,a}' + E_b $$
## Where $T_b$ is an $N \times A$ matrix of block scores, $P_b$ is a
## $K_b \times A$ matrix of block loadings, $E_b$ is an $N \times K_b$ matrix
## of residuals. The number of components $A$ is chosen by cross-validation
## such that the model's cumulative $R^2$ metric is greater than 0.01
## and the model components' cumulative $Q^2$ metrics are strictly
## increasing.
##
## PCA model scores in $T$ (and $T_b$) are low-rank approximations of
## observations in $X$ (and $X_b$), and model loadings in $P$ (and $P_b$)
## are approximations of variables in $X$ (and $X_b$).
## @end tex
## @deftypefn {Function File} {@var{mdl} =} mbpca (@var{X})
## @deftypefnx {Function File} {@var{mdl} =} mbpca (@var{X}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} mbpca (@var{X}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} mbpca (@var{X}, @var{scalefn}, @var{ncv}, @var{aout})
## Performs Multiblock Principal Component Analysis (MBPCA) on the set of
## input data matrices @var{X} = @{@var{X1} @dots{} @var{XB}@} using a
## Consensus PCA (CPCA) method with modifications from Westerhuis. More
## information may be found here:
##
## @quotation
## J. Westerhuis, T. Kourti, J. F. MacGregor. `Analysis of Multiblock and
## Hierarchical PCA and PLS models', Journal of Chemometrics,
## 1998(21): 301-321.
## @end quotation
##
## The optional argument @var{scalefn} may be set to the function handle of
## an appropriate scaling routine. The default scaling function is @ref{suv}.
##
## The optional argument @var{ncv} may be used to specify the number of
## cross-validation iterations and/or groups. If @var{ncv} is a scalar, it
## will be used to set the number of CV groups, with 10 CV iterations. If
## @var{ncv} is a 2-element vector, its elements will be used to set the
## number of CV groups and iterations, respectively. The default is to use
## 7 groups and 10 iterations.
##
## The final optional argument @var{aout} may be used if a specific number of
## model components is desired. @strong{CAUTION:} using @var{aout} will not
## guarantee that the resultant model components are statistically
## significant!
## @end deftypefn

function mdl = mbpca (X, scalefn, ncv, aout)
  % check for proper arguments.
  if (nargin < 1 || !iscell(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a second argument was passed.
  if (nargin < 2 || isempty(scalefn) || !is_function_handle(scalefn))
    % store a default scaling function handle.
    scalefn = @suv;
  end

  % check if a third argument was passed.
  if (nargin < 3 || isempty(ncv))
    % use the default values.
    ncv = [7, 10];
  end

  % check if a target component count was passed.
  if (nargin < 4 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  end

  % ensure the input cell array contains at least two elements.
  if (isempty(X) || length(X) < 2)
    % this won't fly. throw an exception.
    error('mbpca: input data array must contain at least two blocks');
  end

  % ensure the input cell array contains only matrices.
  if (!all(cellfun(@(Xb) ismatrix(Xb), X)))
    % this won't fly either. throw an exception.
    error('mbpca: input data array may contain only matrices');
  end

  % get the number of blocks in the dataset.
  B = length(X);

  % get the number of observations in the dataset.
  N = unique(cellfun(@(Xb) rows(Xb), X));

  % get the number of variables in each block.
  K = reshape(cellfun(@(Xb) columns(Xb), X), B, 1);
  Kmat = [[1; cumsum(K) + 1], [cumsum(K); 1]];

  % ensure only one row count exists in the data.
  if (!isscalar(N))
    % whoops, this really won't fly. throw an exception.
    error('mbpca: blocks must have the same observation count');
  end

  % build the supermatrix: a concatenation of all blocks.
  Xsuper = cell2mat(reshape(X, 1, B));

  % build a superscaling vector that simply normalizes each block's
  % contribution to the variance of the supermatrix.
  wsuper = [];
  for b = 1 : B
    % append the scale factors of the current block.
    wsuper = [wsuper, repmat(sqrt(K(b)), 1, K(b))];
  end

  % build an initial pca model from which to extract block information.
  mdl = pca(Xsuper, scalefn, ncv, aout, wsuper);

  % set the *correct* type of the output structure.
  mdl.type = 'pca';
  mdl.builder = @mbpca;
  mdl.classifier = @mbpcaclassify;

  % initialize the superweights.
  mdl.W = [];

  % initialize the blocks array.
  mdl.B = B;
  mdl.blocks = cell(mdl.B, 1);

  % loop for each block to build the individual block models.
  for b = 1 : B
    % store the block model type.
    mdl.blocks{b}.type = 'pca';
    mdl.blocks{b}.builder = @pca;
    mdl.blocks{b}.classifier = @pcaclassify;

    % store the block scaling information.
    mdl.blocks{b}.scaling = mdl.scaling;

    % store the block data.
    mdl.blocks{b}.X0 = X{b};

    % store the scaled block data.
    [mdl.blocks{b}.X, mdl.blocks{b}.mean.X, mdl.blocks{b}.scale.X] = ...
      mdl.scaling(X{b});

    % store the block size.
    mdl.blocks{b}.N = rows(X{b});
    mdl.blocks{b}.K = columns(X{b});

    % store the block component count.
    mdl.blocks{b}.A = mdl.A;

    % initialize the block component fields.
    mdl.blocks{b}.T = [];
    mdl.blocks{b}.P = [];
    mdl.blocks{b}.E = mdl.blocks{b}.X;

    % initialize more block component fields.
    mdl.blocks{b}.iter = mdl.iter;
    mdl.blocks{b}.lambda = [];
    mdl.blocks{b}.Rsq.X.comp = [];
    mdl.blocks{b}.Qsq.comp = [];

    % loop for each component.
    for a = 1 : mdl.A
      % grab the current component's superscore.
      t = mdl.T(:,a);

      % compute the current component's block loadings.
      pb = (mdl.blocks{b}.E' * t) ./ ((t' * t) * sqrt(K(b)));

      % compute the current component's block scores.
      tb = (mdl.blocks{b}.E * pb) ./ ((pb' * pb) * sqrt(K(b)));

      % compute the block eigenvalue and R-squared value.
      lambda = t' * t;
      Rsq = ssratio(t * pb', mdl.blocks{b}.X ./ sqrt(mdl.blocks{b}.K));

      % compute the block Q-squared value.
      Qsq = 0;

      % store the block eigenvalue, R-squared and Q-squared values.
      mdl.blocks{b}.lambda = [mdl.blocks{b}.lambda; lambda];
      mdl.blocks{b}.Rsq.X.comp = [mdl.blocks{b}.Rsq.X.comp; Rsq];
      mdl.blocks{b}.Qsq.comp = [mdl.blocks{b}.Qsq.comp; Qsq];

      % store the computed scores and loadings.
      mdl.blocks{b}.T = [mdl.blocks{b}.T, tb];
      mdl.blocks{b}.P = [mdl.blocks{b}.P, pb];

      % deflate the block residuals.
      mdl.blocks{b}.E -= t * pb';
    end

    % calculate the block leverages of the observations.
    Tb = mdl.blocks{b}.T;
    mdl.blocks{b}.Hx = diag(Tb * inv(Tb' * Tb) * Tb');

    % store the cumulative block R-squared and Q-squared values.
    mdl.blocks{b}.Rsq.X.cum = sum(mdl.blocks{b}.Rsq.X.comp);
    mdl.blocks{b}.Qsq.cum = sum(mdl.blocks{b}.Qsq.comp);
  end
end

