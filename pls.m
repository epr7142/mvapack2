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
## @anchor{pls}
## @tex
## Partial Least Squares Projection to Latent Structures (PLS) generates
## a linear model of an input data matrix $X$ and response matrix $Y$
## such that the model components are orthogonal and capture a maximal
## amount of correlation between $X$ and $Y$. The PLS model is structured
## as follows:
## $$ X = T P' + E = \sum_{a=1}^A t_a p_a' + E $$
## $$ Y = U C' + F = \sum_{a=1}^A u_a c_a' + F $$
## Where $T$ is an $N \times A$ matrix of scores, $P$ is a $K \times A$
## matrix of $X$-loadings, $E$ is an $N \times K$ matrix of $X$
## residuals, $U$ is an $N \times A$ matrix of $Y$-scores, $C$ is an
## $M \times A$ matrix of $Y$-weights, and $F$ is an $N \times M$
## matrix of $Y$ residuals.
##
## The number of components $A$ is chosen by
## cross-validation such that the model's cumulative $R^2_X$ and
## $R^2_Y$ metrics are greater than 0.01 and the model components'
## cumulative $Q^2$ metrics are strictly increasing.
##
## PLS model scores in $T$ are low-rank approximations of observations
## in $X$ and are good predictors of $U$. Model $X$-loadings in $P$ are
## approximations of variables in $X$ and $Y$-weights in $C$ are
## approximations of columns of $Y$.
## @end tex
## @deftypefn {Function File} {@var{mdl} =} pls (@var{X}, @var{Y})
## @deftypefnx {Function File} {@var{mdl} =} pls (@var{X}, @var{Y}, @var{scalefn})
## @deftypefnx {Function File} {@var{mdl} =} pls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv})
## @deftypefnx {Function File} {@var{mdl} =} pls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout})
## @deftypefnx {Function File} {@var{mdl} =} pls (@var{X}, @var{Y}, @var{scalefn}, @var{ncv}, @var{aout}, @var{w})
## Performs Projection to Latent Structures (PLS) on the input data matrix
## @var{X} and response matrix @var{Y} using the iterative NIPALS method.
## More information may be found here:
##
## @quotation
## S. Wold, M. Sjostrom, L. Eriksson, `PLS-regression: A basic tool in
## chemometrics', Chemometr. Intell. Lab. Syst., 2001(58): 109-130.
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

function mdl = pls (X, Y, scalefn, ncv, aout, w)
  % check for the minimum argument count.
  if (nargin < 2 || !ismatrix(X))
    % print the usage.
    print_usage();
  end

  % set the type of the output structure.
  mdl.type = 'pls';
  mdl.builder = @pls;
  mdl.classifier = @plsclassify;

  % check if the third argument was passed.
  if (nargin >= 3 && !isempty(scalefn) && is_function_handle(scalefn))
    % store the passed scaling function handle.
    mdl.scaling = scalefn;
  else
    % store the default scaling function handle.
    mdl.scaling = @suv;
  end

  % check if a cross-validation argument was passed.
  if (nargin >= 4 && !isempty(ncv))
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
       error('pls: invalid construction of cross-validation argument');
    end
  else
    % no. use the default: 7 groups and 10 iterations.
    mdl.cv.ngroup = 7;
    mdl.cv.niter = 10;
  end

  % check for valid cross-validation parameters.
  if (mdl.cv.ngroup < 2 || mdl.cv.niter < 1)
    % throw an exception.
    error('pls: invalid cross-validation parameters');
  end

  % check if a target component count was passed.
  if (nargin < 5 || isempty(aout))
    % set the target count to zero, indicating a reliance on validation.
    aout = 0;
  end

  % check if a weighting vector was passed.
  if (nargin < 6 || isempty(w))
    % set a unit weight.
    w = ones(1, columns(X));
  end

  % store whether the model is of the discriminant type.
  mdl.isda = isclasses(Y);

  % store the original data and response matrices.
  mdl.X0 = X;
  mdl.Y0 = Y;

  % store the data matrix variable weights.
  mdl.weight.X = w;
  mdl.weight.Y = ones(1, columns(Y));

  % store the scaled data and response matrices.
  [mdl.X, mdl.mean.X, mdl.scale.X] = mdl.scaling(X, w);
  [mdl.Y, mdl.mean.Y, mdl.scale.Y] = mdl.scaling(Y);

  % store the matrix dimensions.
  [mdl.N, mdl.K] = size(mdl.X);
  mdl.M = columns(mdl.Y);
  mdl.A = 0;

  % ensure the cv group count is less than the observation count.
  if (mdl.cv.ngroup > mdl.N)
    % whoops. throw an exception.
    error('pls: number of cross-validation groups exceeds observation count');
  end

  % ensure the observation count matches the response count.
  if (mdl.N != rows(mdl.Y))
    % whoops. throw an exception.
    error('pls: row counts of data and response matrices do not match');
  end

  % initialize the model fields.
  mdl.W = [];
  mdl.T = [];
  mdl.P = [];
  mdl.U = [];
  mdl.C = [];
  mdl.E = mdl.X;
  mdl.F = mdl.Y;

  % initialize eigenvalue, R-squared and Q-squared values.
  mdl.iter = [];
  mdl.lambda = [];
  mdl.Rsq.X.comp = [];
  mdl.Rsq.Y.comp = [];
  mdl.Qsq.comp = [];
  mdl.Qsq.dev = [];

  % initialize the cross-validation structure.
  mdl.CV = cvplsinit(mdl);

  % calculate the maximum possible component count.
  Amax = max([aout; min([mdl.N / 2; mdl.K / 2; mdl.M])]);

  % loop through all possibly significant components.
  for a = 1 : Amax
    % calculate a new component.
    [w, t, p, u, c, iter] = plscomp(mdl.E, mdl.F);

    % calculate the current component contribution.
    TP = t * p';
    TC = t * c';

    % calculate the eigenvalue and R-squared values.
    lambda = t' * t;
    RsqX = ssratio(TP, mdl.X);
    RsqY = ssratio(TC, mdl.Y);

    % calculate the cross-validation Q-squared value.
    [Qsq, Qstd, mdl.CV] = cvpls(mdl);

    % see if we have a target component count.
    if (aout)
      % see if we've reached it.
      if (mdl.A >= aout)
        % we have. break the loop.
        break;
      end
    else
      % check if the current component is significant.
      if (RsqX < 0.01 || RsqY < 0.01 || ...
          sum([mdl.Qsq.comp; Qsq]) < sum(mdl.Qsq.comp))
        % no. break from the loop.
        break;
      end
    end

    % store the new eigenvalue and R-squared values.
    mdl.lambda = [mdl.lambda; lambda];
    mdl.Rsq.X.comp = [mdl.Rsq.X.comp; RsqX];
    mdl.Rsq.Y.comp = [mdl.Rsq.Y.comp; RsqY];
    mdl.Qsq.comp = [mdl.Qsq.comp; Qsq];
    mdl.Qsq.dev = [mdl.Qsq.dev; Qstd];

    % store the new component.
    mdl.W = [mdl.W, w];
    mdl.T = [mdl.T, t];
    mdl.P = [mdl.P, p];
    mdl.U = [mdl.U, u];
    mdl.C = [mdl.C, c];
    mdl.iter = [mdl.iter; iter];

    % store the new un-deflated weights and regression coefficients.
    mdl.Ws = mdl.W * inv(mdl.P' * mdl.W);
    mdl.B = mdl.Ws * mdl.C';

    % store the new residuals and increment the component count.
    mdl.E = mdl.E - TP;
    mdl.F = mdl.F - TC;
    mdl.A++;
  end

  % calculate the leverages of the observations.
  mdl.Hx = diag(mdl.T * inv(mdl.T' * mdl.T) * mdl.T');
  mdl.Hy = diag(mdl.U * inv(mdl.U' * mdl.U) * mdl.U');

  % store the cumulative R-squared and Q-squared values.
  mdl.Rsq.X.cum = sum(mdl.Rsq.X.comp);
  mdl.Rsq.Y.cum = sum(mdl.Rsq.Y.comp);
  mdl.Qsq.cum = sum(mdl.Qsq.comp);

  % check if any significant components were found.
  if (mdl.A == 0)
    % no. throw an exception.
    error('pls: no significant components found');
  end
end

