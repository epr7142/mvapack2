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
## @anchor{mbopls}
## @tex
## Multiblock Orthogonal Projections to Latent Structures (MBOPLS)
## generates a linear model of a set of $B$ input data matrices
## $X = [X_1  X_2 \cdots X_B]$ and a response matrix $Y$ such that the
## overall model components are orthogonal and capture a maximal amount of
## correlation between $X$ and $Y$. The MBOPLS model is structured
## as follows:
## $$ X = [X_1@ @  X_2@ @ \cdots@ @ X_B], \forall b \in [1, B] $$
## Where each block ($X_b$) is itself a bilinear model:
## $$ X_b = T_b P_b' + T_{o,b} P_{o,b}' + E_b = \sum_{a=1}^A t_{b,a} p_{b,a}' + \sum_{a_o=1}^{A_o} t_{o,b,a} p_{o,b,a}' + E_b $$
## $$ Y = U C' + F = \sum_{a=1}^A u_a c_a' + F $$
## Where $T_b$ is an $N \times A$ matrix of block scores, $P_b$ is a
## $K_b \times A$ matrix of $X_b$-loadings, $E_b$ is an $N \times K_b$ matrix
## of $X_b$ residuals, $U$ is an $N \times A$ matrix of $Y$-scores, $C$ is an
## $M \times A$ matrix of $Y$-weights, and $F$ is an $N \times M$
## matrix of $Y$ residuals.
##
## The number of components $A$ is chosen by
## cross-validation such that the model's cumulative $R^2_X$ and
## $R^2_Y$ metrics are greater than 0.01 and the model components'
## cumulative $Q^2$ metrics are strictly increasing.
##
## PLS model scores in $T$ (and $T_b$) are low-rank approximations of
## observations in $X$ (and $X_b$) and are good predictors of $U$. Model
## $X_b$-loadings in $P_b$ are approximations of variables in $X_b$ and
## $Y$-weights in $C$ are approximations of columns of $Y$. Orthogonal
## scores and loadings have been identified in $T_{o,b}$ and $P_{o,b}$
## and are essentially low-rank stores for variation that interferes
## with PLS predictive components.
## @end tex
## @deftypefn {Function File} {@var{mdl} =} mbopls (@var{X}, @var{Y})
## @deftypefnx {Function File} {@var{mdl} =} mbopls (@var{X}, @var{Y}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} mbopls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} mbopls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout})
## Performs Multiblock Orthogonal Projections to Latent Structures (MBOPLS)
## on the set of input data matrices @var{X} = @{@var{X1} @dots{} @var{XB}@}
##and response matrix @var{Y} using a multiblock OPLS method.
## More information may be found here:
##
## @quotation
## B. Worley, R. Powers. `A Sequential Algorithm for Multiblock Orthogonal
## Projections to Latent Structures', Journal of Chemometrics,
## In prep.
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
## model components is desired. If @var{aout} is specified as a scalar, then
## the model will be built with that number of predictive components. If
## @var{aout} is a vector, the first element should hold the desired number of
## predictive components. The remaining elements in @var{aout} should hold the
## desired number of orthogonal components for each predictive component.
## @strong{CAUTION:} using @var{aout} will not guarantee that the resultant
## model components are statistically significant!
## @end deftypefn

function mdl = mbopls (X, Y, scalefn, ncv, aout)
  % check for proper arguments.
  if (nargin < 2 || !iscell(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % check if a second argument was passed.
  if (nargin < 3 || isempty(scalefn) || !is_function_handle(scalefn))
    % store a default scaling function handle.
    scalefn = @spareto;
  end

  % check if a cross-validation argument was passed.
  if (nargin < 4 || isempty(ncv))
    % no. use the default: 7 groups and 10 iterations.
    ncv = [7, 10];
  end

  % check if a target component count was passed.
  if (nargin < 5 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  end

  % ensure the input cell array contains at least two elements.
  if (isempty(X) || length(X) < 2)
    % this won't fly. throw an exception.
    error('mbopls: input data array must contain at least two blocks');
  end

  % ensure the input cell array contains only matrices.
  if (!all(cellfun(@(Xb) ismatrix(Xb), X)))
    % this won't fly either. throw an exception.
    error('mbopls: input data array may contain only matrices');
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
    error('mbopls: blocks must have the same observation count');
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

  % build an initial pls model from which to extract block information.
  mdl = opls(Xsuper, Y, scalefn, ncv, aout, wsuper);

  % set the correct type of the output structure.
  mdl.type = 'opls';
  mdl.builder = @mbopls;
  mdl.classifier = @mboplsclassify;

  % store the original input multiblock array in the model.
  mdl.Xmb = X;

  % rename the regression coefficients so we can use 'B' as a block count.
  mdl.beta = mdl.B;

  % initialize the blocks array.
  mdl.B = B;
  mdl.blocks = cell(mdl.B, 1);

  % loop for each block to build the individual block models.
  for b = 1 : B
    % store the block model type.
    mdl.blocks{b}.type = 'opls';
    mdl.blocks{b}.builder = @opls;
    mdl.blocks{b}.classifier = @oplsclassify;

    % store the block scaling information.
    mdl.blocks{b}.scaling = mdl.scaling;

    % store the block data.
    mdl.blocks{b}.X0 = X{b};
    mdl.blocks{b}.Y0 = mdl.Y0;

    % store the scaled block data.
    [mdl.blocks{b}.X, mdl.blocks{b}.mean.X, mdl.blocks{b}.scale.X] = ...
      mdl.scaling(X{b});

    % store the scaled response data.
    mdl.blocks{b}.Y = mdl.Y;
    mdl.blocks{b}.mean.Y = mdl.mean.Y;
    mdl.blocks{b}.scale.Y = mdl.scale.Y;

    % store the block size.
    mdl.blocks{b}.N = rows(X{b});
    mdl.blocks{b}.K = columns(X{b});

    % store the response count.
    mdl.blocks{b}.M = mdl.M;

    % store the block component count.
    mdl.blocks{b}.A = mdl.A;
    mdl.blocks{b}.Ap = mdl.Ap;
    mdl.blocks{b}.Ao = mdl.Ao;
    mdl.blocks{b}.no = mdl.no;

    % initialize the block component fields.
    mdl.blocks{b}.Wp = [];
    mdl.blocks{b}.Tp = [];
    mdl.blocks{b}.Pp = [];
    mdl.blocks{b}.Wo = [];
    mdl.blocks{b}.To = [];
    mdl.blocks{b}.Po = [];
    mdl.blocks{b}.E = mdl.blocks{b}.X;
    mdl.blocks{b}.F = mdl.Y;

    % initialize more block component fields.
    mdl.blocks{b}.iter = mdl.iter;
    mdl.blocks{b}.lambda = [];
    mdl.blocks{b}.Rsq.Xp.comp = [];
    mdl.blocks{b}.Rsq.Xo.comp = [];
    mdl.blocks{b}.Rsq.Y.comp = [];
    mdl.blocks{b}.Qsq.comp = [];
    mdl.blocks{b}.Qsq.dev = [];

    % initialize the block cross-validation structure.
    mdl.blocks{b}.cv = mdl.cv;
    mdl.blocks{b}.CV = {};

    % get the variable indices of the current block.
    k1 = Kmat(b,1);
    k2 = Kmat(b,2);

    % loop for each component.
    ao = 1;
    for a = 1 : mdl.Ap
      % grab the current predictive superscores.
      t = mdl.Tp(:,a);
      u = mdl.U(:,a);
      c = mdl.C(:,a);

      % loop for each orthogonal component.
      for aoidx = 1 : mdl.no(a)
        % grab the current orthogonal superscores.
        to = mdl.To(:,ao);

        % compute the orthogonal block loadings.
        pbo = (mdl.blocks{b}.E' * to) ./ ((to' * to) * sqrt(K(b)));

        % compute the orthogonal block weights and scores.
        wbo = mdl.Wo(k1:k2,ao);
        tbo = (mdl.blocks{b}.E * wbo) ./ sqrt(K(b));

        % store the computed scores and loadings.
        mdl.blocks{b}.To = [mdl.blocks{b}.To, tbo];
        mdl.blocks{b}.Po = [mdl.blocks{b}.Po, pbo];
        mdl.blocks{b}.Wo = [mdl.blocks{b}.Wo, wbo];

        % compute and store the R-squared value.
        RsqXo = ssratio(tbo * pbo', mdl.blocks{b}.X ./ sqrt(K(b)));
        mdl.blocks{b}.Rsq.Xo.comp = [mdl.blocks{b}.Rsq.Xo.comp; RsqXo];

        % deflate the block residuals.
        mdl.blocks{b}.E -= to * pbo';

        % increment the orthogonal component count.
        ao++;
      end

      % compute the current component's block weights.
      wb = mdl.blocks{b}.E' * u;
      wb = wb ./ sqrt(wb' * wb);

      % compute the current component's block predictive X-scores
      % and X-loadings.
      tb = (mdl.blocks{b}.E * wb) ./ sqrt(K(b));
      pb = (mdl.blocks{b}.E' * t) ./ ((t' * t) * sqrt(K(b)));

      % compute the block eigenvalue and R-squared value.
      lambda = tb' * tb;
      RsqX = ssratio(tb * wb', mdl.blocks{b}.X ./ sqrt(K(b)));
      RsqY = mdl.Rsq.Y.comp(a);

      % compute the block Q-squared value.
      Qsq = mdl.Qsq.comp(a);
      Qdev = mdl.Qsq.dev(a);

      % store the block eigenvalue, R-squared and Q-squared values.
      mdl.blocks{b}.lambda = [mdl.blocks{b}.lambda; lambda];
      mdl.blocks{b}.Rsq.Xp.comp = [mdl.blocks{b}.Rsq.Xp.comp; RsqX];
      mdl.blocks{b}.Rsq.Y.comp = [mdl.blocks{b}.Rsq.Y.comp; RsqY];
      mdl.blocks{b}.Qsq.comp = [mdl.blocks{b}.Qsq.comp; Qsq];
      mdl.blocks{b}.Qsq.dev = [mdl.blocks{b}.Qsq.dev; Qdev];

      % store the computed scores and loadings.
      mdl.blocks{b}.Wp = [mdl.blocks{b}.Wp, wb];
      mdl.blocks{b}.Tp = [mdl.blocks{b}.Tp, tb];
      mdl.blocks{b}.Pp = [mdl.blocks{b}.Pp, pb];

      % deflate the block residuals.
      mdl.blocks{b}.E -= t * pb';
      mdl.blocks{b}.F -= t * c';
    end

    % concatenate the predictive and orthogonal components.
    mdl.blocks{b}.W = [mdl.blocks{b}.Wp, mdl.blocks{b}.Wo];
    mdl.blocks{b}.T = [mdl.blocks{b}.Tp, mdl.blocks{b}.To];
    mdl.blocks{b}.P = [mdl.blocks{b}.Pp, mdl.blocks{b}.Po];

    % calculate the block leverages of the observations.
    Tb = mdl.blocks{b}.Tp;
    Ub = mdl.U;
    mdl.blocks{b}.Hx = diag(Tb * inv(Tb' * Tb) * Tb');
    mdl.blocks{b}.Hy = diag(Ub * inv(Ub' * Ub) * Ub');

    % store the cumulative block R-squared and Q-squared values.
    mdl.blocks{b}.Rsq.Xp.cum = sum(mdl.blocks{b}.Rsq.Xp.comp);
    mdl.blocks{b}.Rsq.Xo.cum = sum(mdl.blocks{b}.Rsq.Xo.comp);
    mdl.blocks{b}.Rsq.Y.cum = sum(mdl.blocks{b}.Rsq.Y.comp);
    mdl.blocks{b}.Qsq.cum = sum(mdl.blocks{b}.Qsq.comp);
  end
end

