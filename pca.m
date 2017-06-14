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
## @anchor{pca}
## @tex
## Principal Component Analysis (PCA) generates a linear model of an input
## data matrix $X$ such that the model components are orthogonal and
## capture a maximal amount of variation present in $X$. The PCA
## model is structured as follows:
## $$ X = T P' + E = \sum_{a=1}^A t_a p_a' + E $$
## Where $T$ is an $N \times A$ matrix of scores, $P$ is a
## $K \times A$ matrix of loadings, $E$ is an $N \times K$ matrix of
## residuals. The number of components $A$ is chosen by cross-validation
## such that the model's cumulative $R^2$ metric is greater than 0.01
## and the model components' cumulative $Q^2$ metrics are strictly
## increasing.
##
## PCA model scores in $T$ are low-rank approximations of observations
## in $X$, and model loadings in $P$ are approximations of variables
## in $X$.
## @end tex
## @deftypefn {Function File} {@var{mdl} =} pca (@var{X})
## @deftypefnx {Function File} {@var{mdl} =} pca (@var{X}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} pca (@var{X}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} pca (@var{X}, @var{scalefn}, @var{ncv}, @var{aout})
## @deftypefnx {Function File} {@var{mdl} =} pca (@var{X}, @var{scalefn}, @var{ncv}, @var{aout}, @var{w})
## Performs Principal Component Analysis (PCA) on the input data matrix
## @var{X} using the iterative NIPALS method. More information may be found
## here:
##
## @quotation
## P. Geladi, B. R. Kowalski. `Partial Least Squares Regression: A Tutorial',
## Analytica Chimica Acta, 1986(185): 1-17.
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
## The optional argument @var{aout} may be used if a specific number of
## model components is desired. @strong{CAUTION:} using @var{aout} will not
## guarantee that the resultant model components are statistically
## significant!
##
## The last optional argument @var{w} may be used to manually weight the
## variables during the PCA analysis. The weights will be effectively
## multiplied into the scale values obtained by @var{scalefn}.
## @end deftypefn

function mdl = pca (X, scalefn, ncv, aout, w)
  % check for proper arguments.
  if (nargin < 1 || !ismatrix(X))
    % improper arguments. print the usage statement.
    print_usage();
  end

  % set the type of the output structure.
  mdl.type = 'pca';
  mdl.builder = @pca;
  mdl.classifier = @pcaclassify;

  % check if a second argument was passed.
  if (nargin >= 2 && !isempty(scalefn) && is_function_handle(scalefn))
    % store the scaling function handle.
    mdl.scaling = scalefn;
  else
    % store a default scaling function handle.
    mdl.scaling = @suv;
  end

  % check if a cross-validation argument was passed.
  if (nargin >= 3 && !isempty(ncv))
    % vector or scalar?
    if (isvector(ncv) && length(ncv) == 2)
      % use the argument as the groups and iterations.
      mdl.cv.ngroup = ncv(1);
      mdl.cv.niter = ncv(2);
    elseif (isscalar(ncv))
      % use the argument as the number of groups.
      mdl.cv.ngroup = ncv;
      mdl.cv.niter = 10;
    else
      % output an error.
      error('pca: invalid construction of cross-validation argument');
    end
  else
    % no. use the default: 7 groups and 10 iterations.
    mdl.cv.ngroup = 7;
    mdl.cv.niter = 10;
  end

  % check for valid cross-validation parameters.
  if (mdl.cv.ngroup < 2 || mdl.cv.niter < 1)
    % throw an exception.
    error('pca: invalid cross-validation parameters');
  end

  % check if a target component count was passed.
  if (nargin < 4 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  end

  % check if a weighting vector was passed.
  if (nargin < 5 || isempty(w))
    % set a unit weight.
    w = ones(1, columns(X));
  end

  % initialize the discrimination status.
  mdl.isda = false;

  % store the original and scaled data matrix.
  mdl.X0 = X;
  mdl.weight.X = w;
  [mdl.X, mdl.mean.X, mdl.scale.X] = mdl.scaling(X, w);

  % store the matrix dimensions.
  [mdl.N, mdl.K] = size(mdl.X);

  % ensure the cv group count is less than the observation count.
  if (mdl.cv.ngroup > min(mdl.N, mdl.K))
    % whoops. throw an exception.
    error('pca: number of cross-validation groups is too high');
  end

  % initialize the model fields.
  mdl.A = 0;
  mdl.T = [];
  mdl.P = [];
  mdl.E = mdl.X;

  % initialize more model fields.
  mdl.iter = [];
  mdl.lambda = [];
  mdl.Rsq.X.comp = [];
  mdl.Qsq.comp = [];
  mdl.Qsq.dev = [];

  % initialize the cross-validation structure.
  mdl.CV = cvpcainit(mdl);

  % calculate the maximum number of components to calculate.
  Amax = max([aout; min([mdl.N / 2; mdl.K / 2])]);

  % loop through all possibly significant components.
  for a = 1 : Amax
    % calculate a new component.
    [t, p, iter] = __pcacomp(mdl.E, rand(mdl.N, 1));

    % calculate the eigenvalue and R-squared values.
    lambda = t' * t;
    Rsq = ssratio(t * p', mdl.X);

    % calculate the cross-validation Q-squared value.
    [Qsq, Qstd, mdl.CV] = cvpca(mdl, t, p);

    % see if we have a target component count.
    if (aout)
      % see if we've reached it.
      if (mdl.A >= aout)
        % we have. break the loop.
        break;
      end
    else
      % see if the newly calculated component is significant.
      if (Rsq < 0.01 || sum([mdl.Qsq.comp; Qsq]) < sum(mdl.Qsq.comp))
        % no. break execution.
        break;
      end
    end

    % store the new component scores and loadings.
    mdl.T = [mdl.T, t];
    mdl.P = [mdl.P, p];
    mdl.E = mdl.E - t * p';

    % store the iteration count.
    mdl.iter = [mdl.iter; iter];

    % store the eigenvalue, R-squared and Q-squared values.
    mdl.lambda = [mdl.lambda; lambda];
    mdl.Rsq.X.comp = [mdl.Rsq.X.comp; Rsq];
    mdl.Qsq.comp = [mdl.Qsq.comp; Qsq];
    mdl.Qsq.dev = [mdl.Qsq.dev; Qstd];

    % increment the component count.
    mdl.A++;
  end

  % check if any components are significant.
  if (mdl.A == 0)
    % no. throw an exception.
    error('pca: no significant components found');
  end

  % calculate the leverages of the observations.
  mdl.Hx = diag(mdl.T * inv(mdl.T' * mdl.T) * mdl.T');

  % store the cumulative R-squared and Q-squared values.
  mdl.Rsq.X.cum = sum(mdl.Rsq.X.comp);
  mdl.Qsq.cum = sum(mdl.Qsq.comp);
end

